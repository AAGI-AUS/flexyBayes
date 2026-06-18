# fb_from_greta -- IR ingest path for direct-greta-model entry.
#
# Counterpart to fb_from_asreml() (R/fb.R) and fb_from_brms()
# (R/fb_from_brms.R). Builds an `fb_terms` intermediate representation
# (IR) post-hoc from a user-supplied greta_model object, recording
# enough structure that the rest of the flexyBayes post-fit surface
# (summary / coef / vcov / prior_summary / triangulate) works on an
# fb_greta() fit.
#
# The flow runs
# in the opposite direction to the formula entries: there is nothing
# to emit (the user has already built the greta graph); the IR is
# built post-hoc by introspecting `model$target_greta_arrays` and
# `model$dag$node_list`.
#
# Internal -- not exported.

# ---------------------------------------------------------------- #
# Introspection helpers                                            #
# ---------------------------------------------------------------- #

# Enumerate the node list of a greta_model. Stable across greta
# 0.5.x. Returns a list of R6 node objects.
.greta_node_list <- function(model) {
  if (is.null(model$dag) || !is.environment(model$dag)) {
    stop(
      "greta model graph not available: `model$dag` is NULL or ",
      "not an environment. Confirm `model` was built via ",
      "`greta::model()`.",
      call. = FALSE
    )
  }
  if (!exists("node_list", envir = model$dag, inherits = FALSE)) {
    stop(
      "greta model graph is missing `node_list`. Likely cause: a ",
      "greta version incompatibility (this code was validated ",
      "against greta 0.5.x).",
      call. = FALSE
    )
  }
  get("node_list", envir = model$dag, inherits = FALSE)
}

# Identify the likelihood node(s): distribution_node instances whose
# target is a data_node (i.e., bound by `distribution(<data>) <- ...`).
# Returns a list with `family` (character; the distribution_name) and
# `n_data` (integer; the data dimension).
.greta_detect_likelihood <- function(model) {
  nl <- .greta_node_list(model)
  is_dist <- vapply(
    nl,
    function(n) inherits(n, "distribution_node"),
    logical(1)
  )

  # A likelihood distribution is one whose `target` field points to a
  # data_node. Variable-node priors also live in distribution_node
  # instances but their target is a variable_node.
  lik_nodes <- list()
  lik_dims <- integer(0)
  for (i in which(is_dist)) {
    tgt <- nl[[i]]$target
    if (!is.null(tgt) && inherits(tgt, "data_node")) {
      lik_nodes <- c(lik_nodes, list(nl[[i]]))
      d <- tgt$dim
      lik_dims <- c(
        lik_dims,
        if (is.null(d)) NA_integer_ else as.integer(d[1L])
      )
    }
  }

  if (length(lik_nodes) == 0L) {
    # No data-bound likelihood. A pure-prior model can still be sampled
    # but has no observed data; we accept it with a structured note.
    return(list(
      family = "none",
      n_data = 0L,
      response_label = "<no observed data>"
    ))
  }
  if (length(lik_nodes) > 1L) {
    stop(
      "Found ",
      length(lik_nodes),
      " data-bound likelihood ",
      "distributions; fb_greta() supports models with a single ",
      "response. Multi-response support is queued for a ",
      "subsequent release.",
      call. = FALSE
    )
  }

  list(
    family = lik_nodes[[1L]]$distribution_name %||% "unknown",
    n_data = lik_dims[1L],
    response_label = "<greta-direct response>"
  )
}

# Coalescing helper (greta's runtime does not export %||%).
`%||%` <- function(x, y) if (is.null(x)) y else x

# Map a greta distribution_name to a base-R family name where the
# semantics agree. flexyBayes downstream methods key off the family
# name; unknown families pass through as character.
.greta_family_to_fb <- function(distname) {
  switch(
    distname,
    normal = "gaussian",
    student = "gaussian", # approximation; user warned downstream
    bernoulli = "binomial",
    binomial = "binomial",
    poisson = "poisson",
    negative_binomial = "negative_binomial",
    gamma = "gamma",
    beta = "beta",
    lognormal = "gaussian", # log-scale; warned downstream
    distname # unknown -> pass through verbatim
  )
}


#' Ingest a user-built greta model into the flexyBayes IR
#'
#' Wraps a `greta_model` (built with `greta::model(...)`) in a `fb_terms`
#' object -- flexyBayes's backend-agnostic intermediate representation
#' (IR) -- so a natively-specified greta graph can flow through the same
#' downstream machinery (summaries, [triangulate()], canonical-name
#' mapping) as a formula-ingested model. Unlike the formula adapters this
#' is a post-hoc wrapper around an already-built graph, so the resulting
#' IR is greta-only by construction: it carries no fixed / random / rcov
#' term lists, only the populated `greta_meta` slot.
#'
#' @param model A `greta_model` returned by `greta::model(...)`.
#' @param data Optional data.frame used to build the model, recorded on
#'   the IR for downstream methods. Not required for fitting -- greta has
#'   already captured the data into its TensorFlow graph.
#' @param prior Optional [fb_prior()] object. When supplied, every name
#'   must resolve to a target greta_array; semantic agreement with the
#'   graph's encoded priors is the caller's responsibility.
#' @param canonical_names Optional named character vector mapping
#'   greta-side parameter names to canonical names. `NULL` falls back to
#'   the verbatim greta names (with a one-time silenceable note, since
#'   [triangulate()] against a non-greta backend needs canonical names).
#' @param known_matrices Named list, mirroring the ASReml entry; recorded
#'   in `data_summary$known_matrices`.
#'
#' @return An `fb_terms` object with `source = "greta"`, `intercept = NA`,
#'   empty fixed / random / rcov term lists, and the populated
#'   `greta_meta` slot (carrying the model graph, target arrays, and the
#'   canonical-name map). Pass the returned IR straight to [fb()] to fit
#'   the graph via greta while keeping a canonical-name map for
#'   [triangulate()]: `fb(fb_from_greta(model, canonical_names = ...))`.
#'
#' @family flexyBayes ingest adapters
#' @seealso [fb()] and [flexybayes()] for the universal fitting entry
#'   (which accept the returned IR directly); [fb_from_asreml()] and
#'   [fb_from_brms()] for the formula dialects.
#' @examples
#' \dontrun{
#' library(greta)
#' mu    <- normal(0, 10)
#' sigma <- normal(0, 5, truncation = c(0, Inf))
#' y     <- as_data(rnorm(20))
#' distribution(y) <- normal(mu, sigma)
#' m  <- model(mu, sigma)
#' ir <- fb_from_greta(m, canonical_names = c(mu = "(Intercept)"))
#' fit <- fb(ir)   # fit the native graph via greta, keeping the map
#' }
#' @export
fb_from_greta <- function(
  model,
  data = NULL,
  prior = NULL,
  canonical_names = NULL,
  known_matrices = list()
) {
  if (!inherits(model, "greta_model")) {
    stop(
      "`model` must be a `greta_model` object built by ",
      "`greta::model(...)`. Got class: ",
      paste(class(model), collapse = "/"),
      ".",
      call. = FALSE
    )
  }

  targets <- model$target_greta_arrays
  if (length(targets) == 0L) {
    stop(
      "`model$target_greta_arrays` is empty: no target parameters ",
      "for inference. Supply at least one greta_array to ",
      "`greta::model(...)`.",
      call. = FALSE
    )
  }

  target_names <- names(targets)
  if (is.null(target_names) || any(!nzchar(target_names))) {
    stop(
      "All target greta_arrays must be named. Supply named ",
      "arguments to `greta::model(...)`.",
      call. = FALSE
    )
  }

  # Canonical-name resolution:
  # supplied map is validated; otherwise falls back to the verbatim
  # greta names with a one-time structured note (silenceable).
  canonical_map <- if (is.null(canonical_names)) {
    if (!isTRUE(getOption("flexyBayes.silence_canonical_names_note", FALSE))) {
      message(
        "Note: no `canonical_names` supplied to fb_greta(); ",
        "using verbatim greta names. `triangulate()` against a ",
        "non-greta backend will require an explicit ",
        "`canonical_names` argument or a registered auto-mapper. ",
        "Silence via ",
        "`options(flexyBayes.silence_canonical_names_note = TRUE)`."
      )
    }
    setNames(target_names, target_names)
  } else {
    if (!is.character(canonical_names) || is.null(names(canonical_names))) {
      stop(
        "`canonical_names` must be a named character vector mapping ",
        "greta-side names to canonical names.",
        call. = FALSE
      )
    }
    missing <- setdiff(names(canonical_names), target_names)
    if (length(missing) > 0L) {
      stop(
        "`canonical_names` references parameter(s) not in ",
        "`model$target_greta_arrays`: ",
        paste(dQuote(missing), collapse = ", "),
        ". Available targets: ",
        paste(dQuote(target_names), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    # Fill in missing targets with verbatim names.
    full <- setNames(target_names, target_names)
    full[names(canonical_names)] <- canonical_names
    full
  }

  # Per-target dimensionality. Sum of element-counts across targets =
  # total free-parameter dimensionality of the model.
  model_dim <- vapply(
    targets,
    function(a) {
      d <- dim(a)
      if (is.null(d)) 1L else as.integer(prod(d))
    },
    integer(1)
  )

  # Likelihood detection -- family + n_data + response_label.
  lik <- .greta_detect_likelihood(model)
  fb_family <- .greta_family_to_fb(lik$family)

  # Prior validation (membership only; semantic match against the
  # graph's encoded priors is enforced at fb_greta() entry, where the
  # user-supplied prior can be diff'd against the variable-node
  # distributions on the target arrays).
  if (!is.null(prior)) {
    if (!inherits(prior, "fb_prior")) {
      stop("`prior` must be an `fb_prior` object (or NULL).", call. = FALSE)
    }
  }

  # Greta-meta slot: arrays + canonical_map +
  # model_dim + n_data. Extras (likelihood family, response label)
  # carried for downstream summary().
  #
  # v0.5.0: the `model` graph itself is recorded on
  # the IR so a greta-source IR is self-contained -- `fb(fb_from_greta(m,
  # canonical_names = ...))` carries everything `.fit_native_greta()`
  # needs (the graph to sample, plus the canonical map) without a
  # second reference to the original model object.
  greta_meta <- list(
    model = model,
    arrays = targets,
    canonical_map = canonical_map,
    model_dim = model_dim,
    n_data = lik$n_data,
    likelihood = lik$family,
    response_label = lik$response_label
  )

  data_summary <- list(
    n = lik$n_data,
    known_matrices = names(known_matrices),
    data_supplied = !is.null(data)
  )

  # `priors` slot: when the user supplies an fb_prior, store it
  # verbatim (so prior_summary() reads it). When NULL, store a
  # sentinel marking the priors as graph-encoded.
  priors <- if (is.null(prior)) {
    list(source = "graph_encoded", legacy = FALSE)
  } else {
    prior
  }

  new_fb_terms(
    response = lik$response_label,
    family = fb_family,
    link = NULL,
    intercept = NA,
    fixed_terms = list(),
    random_terms = list(),
    rcov_terms = list(),
    addition_terms = list(),
    priors = priors,
    data_summary = data_summary,
    capabilities = character(),
    source = "greta",
    greta_meta = greta_meta
  )
}
