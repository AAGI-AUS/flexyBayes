# ---------------------------------------------------------------- #
# The special-form argument contract                               #
# ---------------------------------------------------------------- #
#
# Every ASReml special form in .asreml_function_vocabulary() has a fixed
# argument shape: how many positional arguments it takes, and which named
# arguments it answers to. Until 0.10.0.9000 that shape lived only inside
# each switch() branch of .walk(), expressed as whatever indexing the branch
# happened to do -- and a malformed call fell through whichever way that
# indexing happened to fall.
#
# A grammar fuzz corpus over all 22 forms (129 malformed variants) found four
# distinct ways a malformed call escaped the refusal taxonomy:
#
#   1. an untyped base error      -- fa(), dsum(), str(), pol(a, b), s(a, b, c)
#                                    raised subscriptOutOfBounds or a coercion
#                                    error from the branch's own indexing;
#   2. a silent neighbour         -- dsum(~units), dsum(a), dsum(1) lowered to
#                                    a homogeneous `units` residual, fitting a
#                                    one-variance model where a sectioned one
#                                    was written;
#   3. a silent NA variable       -- fifteen zero-argument forms (us(), vm(),
#                                    ar1(), ...) parsed to an IR node carrying
#                                    var = NA_character_;
#   4. a silently dropped argument -- thirty-two calls carrying more arguments
#                                    than the branch reads had the surplus
#                                    discarded without comment.
#
# Classes 1, 3 and 4 are all the same defect seen from three sides: the branch
# reads the arguments it expects and never asks whether that is what it was
# given. So the contract is declared once, here, and checked once, in .walk(),
# before the switch() dispatches. A branch can then index what its contract
# guarantees.
#
# te() / ti() / t2() are deliberately absent. They refuse at parse with
# tensor_smooth_unsupported whatever their arity, and that refusal names an
# alternative, so it is the more useful of the two.
#
# @noRd
.fb_grammar_contract <- function() {
  one <- function(named = character(0)) {
    list(min_pos = 1L, max_pos = 1L, named = named)
  }
  list(
    # One factor or covariate, and nothing else.
    ar1 = one(),
    ar2 = one(),
    cor = one(),
    corh = one(),
    us = one(),
    id = one(),
    ide = one(),
    lin = one(),
    spl = one(),
    # diag() / idh() take no level: conditioning on selected levels is what
    # at()'s second argument is for, and accepting one here would let a
    # conditioned model wear a diagonal model's name. That intent was a
    # comment in .walk() before it was a contract.
    diag = one(),
    idh = one(),
    # at(f) or at(f, level).
    at = list(min_pos = 1L, max_pos = 2L, named = character(0)),
    # fa(f, k) -- both required; k is validated as a literal count below.
    fa = list(min_pos = 2L, max_pos = 2L, named = character(0)),
    # pol(x) or pol(x, degree).
    pol = list(min_pos = 1L, max_pos = 2L, named = character(0)),
    # dsum(~ units | f) -- one formula body.
    dsum = list(min_pos = 1L, max_pos = 1L, named = character(0)),
    # str(~terms) or str(~terms, ~structure).
    str = list(min_pos = 1L, max_pos = 2L, named = character(0)),
    # The known-covariance carriers take the variable positionally, the
    # matrix positionally or as `V =` or through `cov = fb_cov(...)`, and
    # their keyword carriers by name. The keyword set is read from
    # .STAGE5A_VM_RESERVED_NAMES rather than copied: a second list of the
    # same names is the drift this contract exists to prevent.
    vm = list(
      min_pos = 1L,
      max_pos = 2L,
      named = c("V", .STAGE5A_VM_RESERVED_NAMES)
    ),
    # ped() reads only the dense matrix and the sparse-precision switch --
    # .cov_representation_from_call() is called for ped without chol /
    # precision / blocks, so those names are accepted by vm() alone.
    ped = list(
      min_pos = 1L,
      max_pos = 2L,
      named = c("V", "cov", "use_sparse_precision")
    ),
    # mgcv's s() carries its basis arguments by name (k, bs, by, m, ...), so
    # the named set is open. Exactly one smoothed variable: flexyBayes fits
    # univariate smooths only, and s(x, z) is a multivariate smooth.
    s = list(min_pos = 1L, max_pos = 1L, named = NULL)
  )
}


# .fb_call_arg_split() --- positional and named argument counts of a call.
#
# expr[[1]] is the function head, so the arguments are expr[-1]. A named
# argument is one whose name is non-empty.
#
# @noRd
.fb_call_arg_split <- function(expr) {
  n_args <- length(expr) - 1L
  if (n_args <= 0L) {
    return(list(n_pos = 0L, named = character(0)))
  }
  nms <- names(expr)
  arg_nms <- if (is.null(nms)) rep("", n_args) else nms[-1L]
  arg_nms[is.na(arg_nms)] <- ""
  list(
    n_pos = sum(!nzchar(arg_nms)),
    named = arg_nms[nzchar(arg_nms)]
  )
}


# .fb_check_grammar_arity() --- the contract gate.
#
# Raises a typed refusal when a special form is called with the wrong number
# of positional arguments or a named argument it does not answer to. Returns
# invisible(NULL) when the call satisfies its contract, so .walk() can index
# what the contract guarantees.
#
# @noRd
.fb_check_grammar_arity <- function(fn, expr) {
  contract <- .fb_grammar_contract()[[fn]]
  if (is.null(contract)) {
    return(invisible(NULL))
  }
  split <- .fb_call_arg_split(expr)
  written <- paste(deparse(expr), collapse = "")

  # `fa(site, )` parses in R as a call carrying the empty symbol. Reading it
  # raises "argument is missing, with no default" from wherever it is first
  # touched, which is a base error about our internals rather than a
  # statement about the formula.
  args <- as.list(expr)[-1L]
  if (length(args) && any(vapply(args, .fb_arg_is_empty, logical(1L)))) {
    stop(.fb_refusal_condition(
      reason_code = "grammar_argument_empty",
      message = sprintf(
        paste0(
          "`%s` has an empty argument -- a trailing or doubled comma leaves ",
          "a gap where %s() expects a value. Remove the comma, or write the ",
          "argument it is holding a place for."
        ),
        written,
        fn
      ),
      form = fn,
      written = written
    ))
  }

  # A named argument written twice. The 0.10.0.9000 contract validated the
  # *set* of names and never their multiplicity, which left the same defect
  # alive in two shapes: `vm(g, V = G, V = G)` was accepted and lowered to one
  # vm term, with the surviving copy an accident of how .walk() indexes; and
  # `s(x, k = 5, k = 9)` reached mgcv's own matcher and raised "formal
  # argument 'k' matched by multiple actual arguments" -- a base error about
  # our internals, which is escape class 1 of the four this contract closes.
  # Checked before arity so a doubled name is named as one, and independent of
  # `contract$named` so the open-ended s() set is covered too.
  if (length(split$named)) {
    dup <- unique(split$named[duplicated(split$named)])
    if (length(dup)) {
      stop(.fb_refusal_condition(
        reason_code = "grammar_argument_duplicated",
        message = sprintf(
          paste0(
            "`%s` writes %s twice. %s() reads each named argument once, so ",
            "which copy takes effect is an accident of the parser rather ",
            "than a choice -- write the argument once with the value you ",
            "mean."
          ),
          written,
          paste0("`", dup, "`", collapse = " and "),
          fn
        ),
        form = fn,
        written = written
      ))
    }
  }

  if (split$n_pos < contract$min_pos) {
    # fa() without its rank is the one case worth naming on its own: the
    # missing argument is a number the user has to choose, not a variable
    # they forgot to type.
    if (identical(fn, "fa") && split$n_pos >= 1L) {
      stop(.fb_refusal_condition(
        reason_code = "fa_rank_missing",
        message = sprintf(
          paste0(
            "`%s` does not say how many factors to fit. fa() needs a rank: ",
            "write fa(%s, 1) for a single common factor, fa(%s, 2) for two. ",
            "The rank must be a positive whole number written in the ",
            "formula, not a variable."
          ),
          written,
          .dep(expr, 2),
          .dep(expr, 2)
        ),
        form = fn,
        written = written
      ))
    }
    stop(.fb_refusal_condition(
      reason_code = "grammar_arity_too_few",
      message = sprintf(
        paste0(
          "`%s` is missing an argument. %s() takes %s, and %d %s written. ",
          "A special form with a missing argument used to parse into a term ",
          "with no variable attached, which reaches the engine as a model ",
          "nobody wrote."
        ),
        written,
        fn,
        .fb_arity_phrase(contract),
        split$n_pos,
        if (split$n_pos == 1L) "was" else "were"
      ),
      form = fn,
      written = written
    ))
  }

  if (split$n_pos > contract$max_pos) {
    if (identical(fn, "s")) {
      stop(.fb_refusal_condition(
        reason_code = "multivariate_smooth_unsupported",
        message = sprintf(
          paste0(
            "`%s` smooths more than one variable, and flexyBayes fits ",
            "univariate penalised splines only -- use s(x) or spl(x) per ",
            "smooth dimension. Multivariate smooths are planned for a ",
            "future release."
          ),
          written
        ),
        form = fn,
        written = written
      ))
    }
    stop(.fb_refusal_condition(
      reason_code = "grammar_arity_too_many",
      message = sprintf(
        paste0(
          "`%s` carries more arguments than %s() reads. %s() takes %s, and ",
          "%d were written. The surplus used to be dropped in silence, so a ",
          "structure written into the call never reached the fitted model."
        ),
        written,
        fn,
        fn,
        .fb_arity_phrase(contract),
        split$n_pos
      ),
      form = fn,
      written = written
    ))
  }

  if (!is.null(contract$named) && length(split$named)) {
    unknown <- setdiff(split$named, contract$named)
    if (length(unknown)) {
      stop(.fb_refusal_condition(
        reason_code = "grammar_named_argument_unknown",
        message = sprintf(
          paste0(
            "`%s` names an argument %s() does not answer to: %s. %s ",
            "An unrecognised name used to be ignored, which reads as ",
            "accepted."
          ),
          written,
          fn,
          paste0("`", unknown, "`", collapse = ", "),
          if (length(contract$named)) {
            paste0(
              fn,
              "() answers to ",
              paste0("`", contract$named, "`", collapse = ", "),
              "."
            )
          } else {
            paste0(fn, "() takes positional arguments only.")
          }
        ),
        form = fn,
        written = written
      ))
    }
  }

  invisible(NULL)
}


# .fb_arity_phrase() --- "exactly one argument" / "one or two arguments".
#
# @noRd
.fb_arity_phrase <- function(contract) {
  word <- c("no", "one", "two", "three", "four", "five")
  lo <- contract$min_pos
  hi <- contract$max_pos
  nm <- function(k) {
    if (k + 1L <= length(word)) word[[k + 1L]] else as.character(k)
  }
  if (lo == hi) {
    paste0("exactly ", nm(lo), if (lo == 1L) " argument" else " arguments")
  } else {
    paste0(nm(lo), " or ", nm(hi), " arguments")
  }
}


# .fb_literal_count() --- read a positive whole number written in the formula.
#
# fa()'s rank and pol()'s degree are counts the user writes into the call.
# Both used to be read with as.integer(), which turns TRUE into 1, 1.5 into 1,
# "1" into 1, and raises an untyped coercion error on a symbol. Each of those
# is a different model from the one written, or no model at all.
#
# `symbol` is how the count is spelled in the documentation -- k for fa(),
# degree for pol() -- so the refusal can say "fa() requires k >= 1" in the
# package's own vocabulary. The floor message is the wording this refusal has
# carried since 0.4.0; only the set of inputs that reach it is wider.
#
# Returns the integer. Raises `reason_code` when the argument is not a
# positive whole number literal.
#
# @noRd
.fb_literal_count <- function(expr, i, fn, what, symbol, reason_code) {
  written <- paste(deparse(expr), collapse = "")
  raw <- expr[[i]]

  not_a_count <- function(saw) {
    stop(.fb_refusal_condition(
      reason_code = reason_code,
      message = sprintf(
        paste0(
          "`%s` gives %s as the %s. %s() requires %s to be a positive whole ",
          "number written in the formula -- %s(%s, 1) or %s(%s, 2). A %s ",
          "that is not a literal count used to be coerced into one, which ",
          "fits a different model from the one written."
        ),
        written, saw, what, fn, symbol,
        fn, .dep(expr, 2), fn, .dep(expr, 2), what
      ),
      form = fn,
      written = written
    ))
  }

  below_floor <- function(value) {
    stop(.fb_refusal_condition(
      reason_code = reason_code,
      message = sprintf(
        "%s() requires %s >= 1; got %s = %s.",
        fn, symbol, symbol, value
      ),
      form = fn,
      written = written
    ))
  }

  if (is.name(raw) || is.call(raw)) {
    not_a_count(paste0(
      "`", paste(deparse(raw), collapse = ""), "`, which is not a number,"
    ))
  }
  if (is.character(raw)) {
    not_a_count(paste0("the string \"", raw, "\""))
  }
  if (is.logical(raw)) {
    not_a_count(paste0("`", as.character(raw), "`, a logical value,"))
  }
  if (!is.numeric(raw) || length(raw) != 1L || is.na(raw)) {
    not_a_count("a value that is not a single number")
  }
  if (!is.finite(raw) || raw != round(raw)) {
    not_a_count(paste0("`", format(raw), "`, which is not a whole number,"))
  }
  if (raw < 1L) {
    below_floor(format(raw))
  }
  as.integer(raw)
}


# .fb_arg_is_empty() --- is this argument R's empty symbol?
#
# `f(a, )` parses to a call whose third element is the empty symbol. It is not
# missing in the length() sense and it cannot be touched without raising, so
# it is tested by identity against quote(expr = ).
#
# @noRd
.fb_arg_is_empty <- function(a) {
  identical(a, quote(expr = ))
}
