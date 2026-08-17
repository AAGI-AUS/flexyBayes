# model_fingerprint -- what has to match before two fits can be compared.
#
# triangulate() reports distances between two marginal posteriors. A
# distance is only informative when the two engines were asked the same
# question, so each fit carries a fingerprint of the question: the
# canonical formula triple, the family and link, a content digest of the
# data, and the prior recorded for every variance component (or the record
# that no prior was recorded and the engine chose its own).
#
# The fingerprint is built at fit time, inside emit_brms() and emit_inla(),
# from the IR and the analysis data -- not from the user's call, which two
# grammars spell differently for the same model. It is stored on
# `fit$extras$fingerprint` on both engines.
#
# Nothing here is cryptographic. The data digest is a deterministic
# column-wise summary designed to notice a different dataset, a reordered
# dataset, or a changed column; it is not a defence against a constructed
# collision, and the fingerprint is documented that way.

# ---------------------------------------------------------------- #
# Fingerprint construction                                          #
# ---------------------------------------------------------------- #

# Build the fingerprint for one fit. `prior_vc_sd` is the legacy scalar
# the caller would have used had no fb_prior() been supplied; it is
# recorded so a legacy-bridge fit says which scalar it ran under.
.fb_model_fingerprint <- function(fb, data, prior_vc_sd = NA_real_) {
  forms <- .fb_canonical_model_strings(fb)
  priors <- .fb_prior_record(fb, prior_vc_sd = prior_vc_sd)

  list(
    version = 1L,
    response = fb$response %||% NA_character_,
    fixed = forms$fixed,
    random = forms$random,
    residual = forms$residual,
    family = as.character(fb$family %||% NA_character_)[1L],
    link = as.character(fb$link %||% NA_character_)[1L],
    data_dim = paste0(nrow(data), "x", ncol(data)),
    data_columns = paste(names(data), collapse = ","),
    data_digest = .fb_data_digest(data),
    priors = priors$recorded,
    prior_sources = priors$sources,
    engine_default_params = priors$engine_default
  )
}

# The elements compared, in the order they are reported. `priors` is
# handled separately (per key) so a mismatch names the parameter.
.FB_FINGERPRINT_ELEMENTS <- c(
  "response",
  "fixed",
  "random",
  "residual",
  "family",
  "link",
  "data_dim",
  "data_columns",
  "data_digest"
)

# Human-readable names for the elements above, used in the refusal.
.FB_FINGERPRINT_LABELS <- c(
  response = "the response variable",
  fixed = "the fixed-effect structure",
  random = "the random-effect structure",
  residual = "the residual structure",
  family = "the response family",
  link = "the link function",
  data_dim = "the data dimensions",
  data_columns = "the data column names",
  data_digest = "the data contents"
)

# Compare two fingerprints. Returns NULL when they agree, otherwise a
# list(element, label, a, b) describing the FIRST disagreement in the
# order above -- first because a reader fixes one thing at a time, and
# the earliest difference is usually the cause of the later ones.
.fb_fingerprint_first_difference <- function(fp_a, fp_b) {
  for (el in .FB_FINGERPRINT_ELEMENTS) {
    va <- fp_a[[el]] %||% NA_character_
    vb <- fp_b[[el]] %||% NA_character_
    if (!identical(va, vb)) {
      return(list(
        element = el,
        label = unname(.FB_FINGERPRINT_LABELS[[el]]),
        a = .fb_fingerprint_show(va),
        b = .fb_fingerprint_show(vb)
      ))
    }
  }

  # Priors are keyed by canonical parameter name. Only a parameter both
  # fits recorded can disagree; a parameter one fit never recorded is not
  # a mismatch, it is an exclusion, and the matched-prior gate handles it.
  pa <- fp_a$priors %||% character(0)
  pb <- fp_b$priors %||% character(0)
  shared_keys <- intersect(names(pa), names(pb))
  for (k in sort(shared_keys)) {
    if (!identical(unname(pa[[k]]), unname(pb[[k]]))) {
      return(list(
        element = paste0("prior[", k, "]"),
        label = paste0("the prior on `", k, "`"),
        a = .fb_fingerprint_show(pa[[k]]),
        b = .fb_fingerprint_show(pb[[k]])
      ))
    }
  }

  NULL
}

# Render one fingerprint element for a message. Long digests are
# truncated -- the reader needs to see that two strings differ, not to
# read either of them in full.
.fb_fingerprint_show <- function(x) {
  s <- paste(as.character(x), collapse = ", ")
  if (!nzchar(s)) {
    return("<empty>")
  }
  if (nchar(s) > 60L) {
    return(paste0(substr(s, 1L, 57L), "..."))
  }
  s
}

# ---------------------------------------------------------------- #
# Canonical formula strings                                         #
# ---------------------------------------------------------------- #

# Structural fields only. `label`, `deparse` and the parser's cached
# level vectors are excluded: the first two carry the user's spelling,
# which two grammars write differently for the same model, and the third
# is a data-derived cache that would make the string longer without
# making it more discriminating (the data digest already covers the
# data).
.FB_TERM_FIELD_EXCLUSIONS <- c(
  "label",
  "deparse",
  "expr",
  "call",
  "levels",
  "var_n",
  "n_levels",
  "nlev",
  "data"
)

# One term -> one canonical string, "<type>[field=value, ...]" with the
# fields in alphabetical order so field order in the IR cannot change it.
.fb_canonical_term_string <- function(term) {
  ttype <- term$type %||% "<unknown>"
  keep <- setdiff(names(term), c("type", .FB_TERM_FIELD_EXCLUSIONS))
  keep <- keep[nzchar(keep)]
  parts <- character(0)
  for (nm in sort(keep)) {
    v <- term[[nm]]
    if (is.null(v) || !is.atomic(v) || length(v) == 0L || length(v) > 16L) {
      next
    }
    parts <- c(
      parts,
      paste0(nm, "=", paste(as.character(v), collapse = ","))
    )
  }
  paste0(ttype, "[", paste(parts, collapse = ", "), "]")
}

# The three formula strings. Term order is sorted away: a model is the
# same model whichever order its random effects were written in.
.fb_canonical_model_strings <- function(fb) {
  one_side <- function(terms) {
    if (is.null(terms) || length(terms) == 0L) {
      return("<none>")
    }
    paste(
      sort(vapply(terms, .fb_canonical_term_string, character(1))),
      collapse = " + "
    )
  }
  # An absent residual side and a bare `units` residual are the same
  # model: the ASReml ingest writes the default residual out and the
  # brms-style ingest leaves it implicit. Normalising them together is
  # what keeps the fingerprint a property of the model rather than of the
  # grammar it was written in.
  resid <- fb$residual_terms %||% list()
  is_plain_units <- length(resid) == 1L &&
    identical(resid[[1L]]$type %||% "", "units") &&
    length(setdiff(names(resid[[1L]]), "type")) == 0L
  residual <- if (length(resid) == 0L || is_plain_units) {
    "units"
  } else {
    one_side(resid)
  }

  list(
    fixed = paste0(
      if (isTRUE(fb$intercept)) "1 + " else "0 + ",
      one_side(fb$fixed_terms)
    ),
    random = one_side(fb$random_terms),
    residual = residual
  )
}

# ---------------------------------------------------------------- #
# Data digest                                                       #
# ---------------------------------------------------------------- #

# Column-wise, order-sensitive content digest. Each column contributes a
# handful of vectorised summaries; the order-sensitive term (sum of the
# values weighted by row position) is what makes a permuted dataset read
# as a different dataset, which matters because a random effect indexed
# by row position is a different model after a shuffle.
.fb_data_digest <- function(data) {
  if (!is.data.frame(data) || ncol(data) == 0L) {
    return(NA_character_)
  }
  per_col <- vapply(
    seq_along(data),
    function(j) .fb_column_digest(names(data)[[j]], data[[j]]),
    character(1)
  )
  paste(per_col, collapse = "|")
}

.fb_column_digest <- function(nm, x) {
  n <- length(x)
  n_na <- sum(is.na(x))
  cls <- paste(class(x), collapse = "/")
  w <- seq_len(n)

  body <- if (is.factor(x)) {
    v <- as.integer(x)
    v[is.na(v)] <- 0L
    lv <- levels(x)
    sprintf(
      "nlev=%d;lvchars=%d;first=%s;last=%s;%s",
      length(lv),
      sum(nchar(lv)),
      if (length(lv) > 0L) lv[[1L]] else "",
      if (length(lv) > 0L) lv[[length(lv)]] else "",
      .fb_numeric_summary(as.numeric(v), w)
    )
  } else if (is.character(x)) {
    nc <- nchar(x)
    nc[is.na(nc)] <- 0L
    sprintf("chars=%s", .fb_numeric_summary(as.numeric(nc), w))
  } else if (is.numeric(x) || is.logical(x)) {
    v <- as.numeric(x)
    v[is.na(v)] <- 0
    .fb_numeric_summary(v, w)
  } else {
    v <- suppressWarnings(as.numeric(x))
    if (all(is.na(v))) {
      "opaque"
    } else {
      v[is.na(v)] <- 0
      .fb_numeric_summary(v, w)
    }
  }

  sprintf("%s:%s:n=%d:na=%d:%s", nm, cls, n, n_na, body)
}

# Five numbers per column: total, total of squares, range, and the
# position-weighted total. Formatted at twelve significant digits, which
# is well inside double precision for identical inputs and well outside
# the range where two different columns coincide by accident.
.fb_numeric_summary <- function(v, w) {
  if (length(v) == 0L) {
    return("empty")
  }
  sprintf(
    "s=%.12g;q=%.12g;mn=%.12g;mx=%.12g;p=%.12g",
    sum(v),
    sum(v * v),
    min(v),
    max(v),
    sum(v * w)
  )
}

# ---------------------------------------------------------------- #
# Prior record                                                      #
# ---------------------------------------------------------------- #

# What prior each variance component actually carries, keyed by the
# canonical parameter name triangulate() compares on.
#
# Three sources, and the difference between them is the whole point:
#
#   "shared_default"  the package's own uniform-on-SD default, injected
#                     identically for every engine before dispatch
#   "explicit"        a user-supplied fb_prior() spec
#   "engine_default"  no prior this package chose -- the term is outside
#                     the default-prior walker, or the run used the legacy
#                     scalar bridge, which reaches brms and not INLA
#
# Only the first two are comparable across engines. The third is recorded
# rather than dropped so that the gate can say *why* a parameter was left
# out instead of quietly comparing two different questions.
.fb_prior_record <- function(fb, prior_vc_sd = NA_real_) {
  targets <- .fb_default_prior_targets(fb)
  engine_default <- stats::setNames(
    targets$engine$reason,
    targets$engine$parameter
  )

  prior <- fb$priors
  if (!inherits(prior, "fb_prior")) {
    # Legacy scalar bridge (or no prior at all). brms receives
    # lognormal(0, prior_vc_sd) on every variance component and INLA
    # receives nothing, keeping its own loggamma default, so no variance
    # component is matched across engines under this route.
    unmatched <- c(
      "sigma",
      paste0("sd_", c(targets$shared, targets$vm_ped))
    )
    reason <- if (is.finite(prior_vc_sd)) {
      paste0(
        "no fb_prior() was supplied: the brms path used the legacy scalar ",
        "lognormal(0, ", format(prior_vc_sd, digits = 6), ") and the INLA ",
        "path kept its own default hyperprior"
      )
    } else {
      paste0(
        "no fb_prior() was supplied, so each engine used its own default ",
        "hyperprior"
      )
    }
    engine_default <- c(
      stats::setNames(rep(reason, length(unmatched)), unmatched),
      engine_default
    )
    return(list(
      recorded = character(0),
      sources = character(0),
      engine_default = engine_default[!duplicated(names(engine_default))]
    ))
  }

  is_default <- !is.null(attr(prior, "fb_prior_default_basis"))
  src <- if (is_default) "shared_default" else "explicit"

  recorded <- character(0)
  sources <- character(0)
  for (s in prior$specs) {
    key <- .fb_prior_spec_parameter(s)
    if (is.na(key)) {
      next
    }
    recorded[[key]] <- .fb_prior_spec_string(s)
    sources[[key]] <- src
  }

  # A term the walker skipped may still have been given an explicit prior
  # by hand; the recorded entry wins over the engine-default note.
  engine_default <- engine_default[
    !(names(engine_default) %in% names(recorded))
  ]

  list(
    recorded = recorded,
    sources = sources,
    engine_default = engine_default
  )
}

# fb_prior target -> canonical parameter name. Fixed-effect and smooth
# targets return NA: they are not variance components, and the fixed part
# of the model is already compared through the formula strings.
.fb_prior_spec_parameter <- function(s) {
  ttype <- s$target$type %||% ""
  if (identical(ttype, "sigma")) {
    return("sigma")
  }
  if (identical(ttype, "sd")) {
    grp <- s$target$group %||% ""
    if (!nzchar(grp)) {
      return(NA_character_)
    }
    return(paste0("sd_", grp))
  }
  NA_character_
}

# Density expression as one canonical string, arguments in alphabetical
# order so two spellings of the same prior compare equal.
.fb_prior_spec_string <- function(s) {
  fam <- s$spec$family %||% "<unknown>"
  args <- s$spec$args %||% list()
  if (length(args) == 0L) {
    return(paste0(fam, "()"))
  }
  nms <- sort(names(args))
  parts <- vapply(
    nms,
    function(a) {
      v <- args[[a]]
      paste0(
        a,
        "=",
        if (is.numeric(v)) sprintf("%.12g", v) else as.character(v)
      )
    },
    character(1)
  )
  paste0(fam, "(", paste(parts, collapse = ", "), ")")
}

# ---------------------------------------------------------------- #
# Accessors                                                         #
# ---------------------------------------------------------------- #

# The fingerprint of a fit, or NULL when the object carries none --
# a fit imported from another package, an aggregated-likelihood fit, or
# anything a user built by hand for a triangulate() peer method.
.fb_fit_fingerprint <- function(fit) {
  if (!is.list(fit)) {
    return(NULL)
  }
  fit$extras$fingerprint
}
