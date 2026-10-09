# model_manifest.R -- the compiled model manifest.
#
# Every fit records, at `fit$extras$manifest`, which stated or defaulted prior
# became which engine slot, the prior the engine reports it received, the
# route the fit ran on, and whether its marginal likelihood is available.
# Reports of a prior, a route or an evidence value read this record.
#
#   inputs    one row per prior specification, stated or defaulted: target,
#             source, and the number of engine slots that consumed it
#   slots     one row per engine slot: what the package sent, the input it
#             came from, and what the engine reports after fitting
#   route     the representation the fit ran on
#   evidence  the marginal likelihood, or the reason it is not available
#
# The read-back is the engine's own record: `result$all.hyper` on INLA, which
# INLA::inla.priors.used() prints, and brms::prior_summary() on brms. A slot
# the package sent no prior to is compared with the default the engine
# declares for it (INLA::inla.models(), INLA::inla.set.control.fixed.default()
# and the `source` column brms writes); a slot matching neither is recorded as
# unexpected.

# The sources an input can have, in the vocabulary .fb_prior_record()
# (R/model_fingerprint.R) records.
.FB_MANIFEST_INPUT_SOURCES <- c(
  "explicit",
  "shared_default",
  "legacy_scalar",
  "filled_default",
  "covariate_scaled_default",
  "structure_default",
  "extra_prior"
)

# The package-chosen sources, as against a prior the caller wrote.
.FB_MANIFEST_DEFAULT_SOURCES <- c(
  "shared_default",
  "filled_default",
  "covariate_scaled_default",
  "structure_default"
)


# --- accessors ------------------------------------------------------------- #

#' The model manifest of a fit
#'
#' @param fit A flexyBayes fit.
#' @returns The `fb_manifest` list, or `NULL` for an object that carries none.
#' @noRd
#' @keywords internal
.fb_manifest <- function(fit) {
  if (!is.list(fit) || !is.list(fit[["extras"]])) {
    return(NULL)
  }
  fit[["extras"]][["manifest"]]
}

#' Where a resolved `fb_prior` came from
#'
#' @param prior An `fb_prior` object.
#' @returns `"shared_default"` for the package's auto-default, which carries
#'   the `fb_prior_default_basis` attribute, and `"explicit"` otherwise.
#' @noRd
#' @keywords internal
.fb_manifest_prior_source <- function(prior) {
  if (is.null(attr(prior, "fb_prior_default_basis"))) {
    "explicit"
  } else {
    "shared_default"
  }
}

#' Whether a manifest's prior was the caller's or the package's
#'
#' @param manifest An `fb_manifest`.
#' @returns `"user"` when any input was written by the caller, `"auto"` when
#'   every input is a package default, and `NA` otherwise.
#' @noRd
#' @keywords internal
.fb_manifest_origin <- function(manifest) {
  sources <- manifest[["inputs"]][["source"]]
  if (any(sources == "explicit")) {
    return("user")
  }
  if (any(sources %in% .FB_MANIFEST_DEFAULT_SOURCES)) {
    return("auto")
  }
  NA_character_
}


# --- inputs ---------------------------------------------------------------- #

.fb_manifest_input_frame <- function(
  target,
  target_key,
  distribution,
  parameters,
  source
) {
  data.frame(
    input_id = paste0(target, ":", target_key, recycle0 = TRUE),
    target = target,
    target_key = target_key,
    distribution = distribution,
    parameters = parameters,
    source = source,
    stringsAsFactors = FALSE
  )
}

#' The key a prior specification is written under
#'
#' @param spec One element of an `fb_prior`'s `specs`.
#' @returns The residual (`"sigma"`), a group, a coefficient or a smoothed
#'   variable, as the caller or the default wrote it.
#' @noRd
#' @keywords internal
.fb_manifest_spec_key <- function(spec) {
  target <- spec[["target"]]
  key <- switch(
    target[["type"]] %||% "",
    sigma = "sigma",
    sd = target[["group"]],
    cor = target[["group"]],
    b = target[["name"]],
    name = target[["name"]],
    smooth = target[["var"]],
    loading = target[["group"]],
    NULL
  )
  if (is.null(key) || !nzchar(key)) NA_character_ else as.character(key)
}

#' The input identifier of a prior specification
#'
#' @param spec One element of an `fb_prior`'s `specs`.
#' @returns `"<target>:<key>"`, the identifier the manifest's inputs carry.
#' @noRd
#' @keywords internal
.fb_manifest_spec_input_id <- function(spec) {
  paste0(spec[["target"]][["type"]], ":", .fb_manifest_spec_key(spec))
}

.fb_manifest_spec_frame <- function(spec, source) {
  text <- .fb_prior_spec_string(spec)
  .fb_manifest_input_frame(
    target = spec[["target"]][["type"]],
    target_key = .fb_manifest_spec_key(spec),
    distribution = spec[["spec"]][["family"]] %||% NA_character_,
    parameters = sub("^[^(]*\\((.*)\\)$", "\\1", text),
    source = source
  )
}

#' The semantic layer: every prior specification a fit was given
#'
#' @param fb An `fb_terms` object.
#' @param scalars A named list of the scalar prior arguments the route used,
#'   `prior_fixed_sd` and `prior_vc_sd`, each a prior density string.
#' @param completion The INLA completion record (R/inla_default_priors.R),
#'   or `NULL`.
#' @param fixed_default_sd The covariate-scaled fixed-effect default SDs,
#'   keyed by coefficient.
#' @param extra A data frame of further inputs in the same columns, or
#'   `NULL`.
#' @returns A data frame with one row per input.
#' @noRd
#' @keywords internal
.fb_manifest_inputs <- function(
  fb,
  scalars = list(),
  completion = NULL,
  fixed_default_sd = numeric(0),
  extra = NULL
) {
  rows <- list(.fb_manifest_input_frame(
    character(0), character(0), character(0), character(0), character(0)
  ))
  prior <- fb[["priors"]]

  # Stated and auto-default specifications ---------------------------------
  if (inherits(prior, "fb_prior")) {
    source <- .fb_manifest_prior_source(prior)
    for (spec in prior[["specs"]]) {
      rows[[length(rows) + 1L]] <- .fb_manifest_spec_frame(spec, source)
    }
  }

  # Scalar prior arguments -------------------------------------------------
  # `prior_vc_sd` is an input whenever the caller wrote it, including beside
  # an fb_prior(), where no route reads it.
  if (
    .fb_prior_scalar_supplied(fb, "vc_sd") &&
      !"prior_vc_sd" %in% names(scalars)
  ) {
    scalars[["prior_vc_sd"]] <- sprintf(
      "lognormal(0, %.12g)",
      .fb_prior_scalar_value(fb, "vc_sd")
    )
  }
  for (nm in names(scalars)) {
    density <- scalars[[nm]]
    rows[[length(rows) + 1L]] <- .fb_manifest_input_frame(
      target = "scalar",
      target_key = nm,
      distribution = sub("\\(.*$", "", density),
      parameters = sub("^[^(]*\\((.*)\\)$", "\\1", density),
      source = "legacy_scalar"
    )
  }

  # What an INLA fit added to the statement --------------------------------
  for (spec in completion[["filled"]] %||% list()) {
    rows[[length(rows) + 1L]] <- .fb_manifest_spec_frame(spec, "filled_default")
  }
  if (!is.null(extra)) {
    rows[[length(rows) + 1L]] <- extra
  }
  inputs <- do.call(rbind, rows)

  # A default bound divided by its covariate's SD replaces the row it came
  # from, because the divided bound is the one the engine received.
  for (spec in completion[["covariate_scaled"]] %||% list()) {
    scaled <- .fb_manifest_spec_frame(spec, "covariate_scaled_default")
    hit <- which(inputs$input_id == scaled$input_id)
    if (length(hit)) {
      inputs[hit[[1L]], ] <- scaled
    } else {
      inputs <- rbind(inputs, scaled)
    }
  }
  for (nm in names(fixed_default_sd)) {
    inputs <- rbind(inputs, .fb_manifest_input_frame(
      target = "b",
      target_key = nm,
      distribution = "normal",
      parameters = sprintf("mean=0, sd=%.12g", fixed_default_sd[[nm]]),
      source = "covariate_scaled_default"
    ))
  }

  inputs$input_id <- make.unique(inputs$input_id, sep = "#")
  rownames(inputs) <- NULL
  inputs
}


# --- slots ----------------------------------------------------------------- #

.fb_manifest_param_string <- function(param) {
  if (!length(param)) {
    return(NA_character_)
  }
  paste(sprintf("%.17g", as.numeric(param)), collapse = ", ")
}

.fb_manifest_param_values <- function(text) {
  if (is.na(text) || !nzchar(text)) {
    return(numeric(0))
  }
  as.numeric(strsplit(text, ", ", fixed = TRUE)[[1L]])
}

# A prior as one line of text: an INLA expression prior is its own text, a
# fixed effect's normal prior is INLA's mean and precision, and any other
# named prior is followed by its parameters.
.fb_manifest_prior_text <- function(
  prior,
  param_text,
  section = NA_character_
) {
  if (is.na(prior) && is.na(param_text)) {
    return(NA_character_)
  }
  if (!is.na(prior) && startsWith(prior, "expression:")) {
    return(prior)
  }
  values <- vapply(
    .fb_manifest_param_values(param_text),
    function(v) format(signif(v, 8L), trim = TRUE),
    character(1L)
  )
  if (identical(section, "fixed") && length(values) == 2L) {
    return(paste0(
      "normal(mean = ", values[[1L]], ", precision = ", values[[2L]], ")"
    ))
  }
  paste0(
    if (is.na(prior)) "" else prior,
    if (length(values)) paste0("(", paste(values, collapse = ", "), ")") else ""
  )
}

#' One engine slot the INLA formula builder wrote
#'
#' @param tag The `f()` term's index name, which INLA reports as its tag.
#' @param model The INLA latent model.
#' @param key The prior key the term reads, or `NA`.
#' @param structure The random-term type.
#' @param entry The hyper entry sent, or `NULL` for none.
#' @param group_model The `control.group` model, or `NA`.
#' @param pinned `TRUE` when the entry writes out the engine's own default.
#' @returns A one-row data frame.
#' @noRd
#' @keywords internal
.fb_manifest_inla_slot <- function(
  tag,
  model,
  key = NA_character_,
  structure = NA_character_,
  entry = NULL,
  group_model = NA_character_,
  pinned = FALSE,
  section = "random"
) {
  data.frame(
    section = section,
    tag = tag,
    model = model,
    target_key = key %||% NA_character_,
    structure = structure,
    emitted_prior = entry[["prior"]] %||% NA_character_,
    emitted_param = .fb_manifest_param_string(entry[["param"]]),
    group_model = group_model,
    pinned = pinned,
    stringsAsFactors = FALSE
  )
}

.fb_manifest_bind_slots <- function(rows) {
  template <- data.frame(
    section = character(0),
    tag = character(0),
    model = character(0),
    target_key = character(0),
    structure = character(0),
    emitted_prior = character(0),
    emitted_param = character(0),
    group_model = character(0),
    pinned = logical(0),
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, c(list(template), rows))
  rownames(out) <- NULL
  out
}

.fb_manifest_bind_reads <- function(rows) {
  template <- data.frame(
    key = character(0),
    structure = character(0),
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, c(list(template), rows))
  rownames(out) <- NULL
  out
}

#' The scalar prior arguments an INLA fit read
#'
#' @param fb An `fb_terms` object.
#' @returns A named list of prior densities: `prior_fixed_sd` when supplied,
#'   and `prior_vc_sd` when supplied without an `fb_prior()`.
#' @noRd
#' @keywords internal
.fb_manifest_inla_scalars <- function(fb) {
  out <- list()
  if (.fb_prior_scalar_supplied(fb, "fixed_sd")) {
    out[["prior_fixed_sd"]] <- sprintf(
      "normal(0, %.12g)",
      .fb_prior_scalar_value(fb, "fixed_sd")
    )
  }
  if (
    .fb_prior_scalar_supplied(fb, "vc_sd") &&
      !inherits(fb[["priors"]], "fb_prior")
  ) {
    out[["prior_vc_sd"]] <- sprintf(
      "lognormal(0, %.12g)",
      .fb_prior_scalar_value(fb, "vc_sd")
    )
  }
  out
}

#' The likelihood slots of an INLA fit
#'
#' @param hyper_ctrl The hyper list as sent.
#' @param inla_family The family, one entry per likelihood.
#' @returns A data frame with one row per likelihood; `target_key` is
#'   `"sigma"` where the residual prior reached it.
#' @noRd
#' @keywords internal
.fb_manifest_inla_family_slots <- function(hyper_ctrl, inla_family) {
  families <- as.character(inla_family)
  entry <- hyper_ctrl[["sigma"]]
  keyword <- .fb_inla_hyper_keyword(families[[1L]])
  carried <- !is.null(entry) && !is.null(keyword)
  rows <- lapply(seq_along(families), function(i) {
    .fb_manifest_inla_slot(
      tag = paste0("INLA.Data", i),
      model = families[[i]],
      key = if (carried) "sigma" else NA_character_,
      structure = if (length(families) > 1L) {
        "sectioned_residual"
      } else {
        "residual"
      },
      entry = if (carried) entry else NULL,
      section = "family"
    )
  })
  .fb_manifest_bind_slots(rows)
}

#' The fixed-effect slots of an INLA fit, read against INLA's own labels
#'
#' @param labels The fixed-effect labels INLA reports.
#' @param control_fixed The `control.fixed` list as sent.
#' @returns A data frame with one row per coefficient that received a prior.
#' @noRd
#' @keywords internal
.fb_manifest_inla_fixed_slots <- function(labels, control_fixed) {
  prec <- control_fixed[["prec"]]
  mean <- control_fixed[["mean"]]
  at <- function(v, nm) {
    if (is.list(v)) v[[nm]] else NULL
  }
  rows <- list()
  for (label in labels) {
    if (identical(label, "(Intercept)")) {
      if (is.null(control_fixed[["prec.intercept"]])) {
        next
      }
      key <- "(Intercept)"
      sent <- c(
        control_fixed[["mean.intercept"]] %||% 0,
        control_fixed[["prec.intercept"]]
      )
    } else if (!is.null(at(prec, label))) {
      key <- label
      sent <- c(at(mean, label) %||% 0, at(prec, label))
    } else if (!is.null(at(prec, "default"))) {
      key <- "default"
      sent <- c(at(mean, "default") %||% 0, at(prec, "default"))
    } else if (is.numeric(prec) && length(prec) == 1L) {
      key <- "default"
      sent <- c(if (is.numeric(mean)) mean else 0, prec)
    } else {
      next
    }
    rows[[length(rows) + 1L]] <- .fb_manifest_inla_slot(
      tag = label,
      model = "fixed",
      key = key,
      structure = "fixed_effect",
      entry = list(prior = "normal", param = sent),
      section = "fixed"
    )
  }
  .fb_manifest_bind_slots(rows)
}


# --- INLA read-back -------------------------------------------------------- #

.fb_manifest_readback_rows <- function(section, tag, hyper) {
  rows <- list()
  for (theta in names(hyper)) {
    h <- hyper[[theta]]
    if (isTRUE(h[["fixed"]])) {
      next
    }
    rows[[length(rows) + 1L]] <- data.frame(
      section = section,
      tag = tag,
      theta = theta,
      short_name = h[["short.name"]] %||% NA_character_,
      prior = as.character(h[["prior"]] %||% NA_character_)[[1L]],
      param = .fb_manifest_param_string(h[["param"]]),
      stringsAsFactors = FALSE
    )
  }
  rows
}

#' The priors an INLA fit reports it used
#'
#' Reads `result$all.hyper` under INLA::inla.priors.used()'s own rule: every
#' section but `predictor`, and only hyperparameters that are not fixed. A
#' `group`, `link` or `mix` block is kept only when INLA's mode names a
#' hyperparameter for it, because `all.hyper` lists one for every term.
#'
#' @param inla_fit An `inla` object.
#' @returns A list: `rows`, a data frame of section, tag, theta, prior and
#'   parameters, with the fixed effects last; and `complete`, whether the
#'   hyperparameter rows number exactly the hyperparameters in INLA's mode.
#' @noRd
#' @keywords internal
.fb_manifest_inla_readback <- function(inla_fit) {
  hyper <- inla_fit[["all.hyper"]]
  theta_names <- names(inla_fit[["mode"]][["theta"]]) %||% character(0)
  rows <- list()
  for (section in c("family", "random")) {
    comps <- hyper[[section]] %||% list()
    for (i in seq_along(comps)) {
      comp <- comps[[i]]
      tag <- as.character(unlist(comp[["hyperid"]] %||% comp[["label"]]))[[1L]]
      rows <- c(rows, .fb_manifest_readback_rows(section, tag, comp[["hyper"]]))
      extra <- list(
        group = comp[["group.hyper"]],
        link = comp[["link"]][["hyper"]],
        mix = comp[["mix"]][["hyper"]]
      )
      prefix <- c(group = "Group", link = "Link", mix = "Mix")
      for (nm in names(extra)) {
        used <- any(startsWith(theta_names, prefix[[nm]]) &
          (section == "family" | endsWith(theta_names, paste0(" for ", tag))))
        if (length(extra[[nm]]) && used) {
          rows <- c(rows, .fb_manifest_readback_rows(nm, tag, extra[[nm]]))
        }
      }
    }
  }
  n_hyper <- length(rows)
  for (fx in hyper[["fixed"]] %||% list()) {
    rows[[length(rows) + 1L]] <- data.frame(
      section = "fixed",
      tag = as.character(fx[["label"]]),
      theta = "beta",
      short_name = NA_character_,
      prior = "normal",
      param = .fb_manifest_param_string(
        c(fx[["prior.mean"]], fx[["prior.prec"]])
      ),
      stringsAsFactors = FALSE
    )
  }
  template <- data.frame(
    section = character(0), tag = character(0), theta = character(0),
    short_name = character(0), prior = character(0), param = character(0),
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, c(list(template), rows))
  rownames(out) <- NULL
  list(rows = out, complete = n_hyper == length(theta_names))
}

#' The prior INLA declares for a slot nobody wrote one for
#'
#' @param section `"family"`, `"random"`, `"group"` or `"fixed"`.
#' @param model The likelihood, latent or group model; for `"fixed"`, the
#'   coefficient label.
#' @param theta The hyperparameter's name within its model.
#' @returns A list with `prior` and `param`, or `NULL` when INLA declares
#'   nothing for it.
#' @noRd
#' @keywords internal
.fb_manifest_inla_default <- function(section, model, theta) {
  if (identical(section, "fixed")) {
    defaults <- INLA::inla.set.control.fixed.default()
    if (identical(model, "(Intercept)")) {
      return(list(
        prior = "normal",
        param = c(defaults[["mean.intercept"]], defaults[["prec.intercept"]])
      ))
    }
    return(list(
      prior = "normal",
      param = c(defaults[["mean"]], defaults[["prec"]])
    ))
  }
  models <- INLA::inla.models()
  block <- switch(
    section,
    family = models[["likelihood"]][[model]][["hyper"]],
    random = models[["latent"]][[model]][["hyper"]],
    group = models[["group"]][[model]][["hyper"]],
    NULL
  )
  entry <- block[[theta]]
  if (is.null(entry)) {
    return(NULL)
  }
  list(prior = as.character(entry[["prior"]])[[1L]], param = entry[["param"]])
}

.fb_manifest_same_prior <- function(prior_a, param_a, prior_b, param_b) {
  if (!is.na(prior_a) && !identical(prior_a, prior_b)) {
    return(FALSE)
  }
  a <- .fb_manifest_param_values(param_a)
  b <- .fb_manifest_param_values(param_b)
  if (!length(a)) {
    return(TRUE)
  }
  length(a) == length(b) && isTRUE(all.equal(a, b, tolerance = 1e-9))
}

# The input a slot's key came from: the stated or defaulted specification
# written under that key, or the scalar argument that covers it.
.fb_manifest_input_for_key <- function(section, key, inputs) {
  if (is.na(key)) {
    return(NA_character_)
  }
  ids <- inputs$input_id
  if (startsWith(key, "wishart:") || startsWith(key, "loading:")) {
    return(if (key %in% ids) key else NA_character_)
  }
  candidates <- switch(
    section,
    fixed = if (identical(key, "default")) {
      "scalar:prior_fixed_sd"
    } else if (identical(key, "(Intercept)")) {
      c("b:(Intercept)", "b:Intercept", "scalar:prior_fixed_sd")
    } else {
      paste0("b:", key)
    },
    family = c("sigma:sigma", "scalar:prior_vc_sd"),
    c(paste0("sd:", key), paste0("smooth:", key), "scalar:prior_vc_sd")
  )
  hit <- candidates[candidates %in% ids]
  if (length(hit)) hit[[1L]] else NA_character_
}

#' The manifest of an INLA fit
#'
#' @param fb The `fb_terms` object the fit was built from.
#' @param inla_fit The `inla` object.
#' @param formula_slots,family_slots Slots from the formula and likelihood
#'   builders.
#' @param reads The prior keys the formula builder consulted, with the
#'   structure consulting each.
#' @param control_fixed The `control.fixed` list as sent.
#' @param scalars,completion,fixed_default_sd Passed to
#'   `.fb_manifest_inputs()`.
#' @param route The route the fit took.
#' @returns An `fb_manifest`.
#' @noRd
#' @keywords internal
.fb_manifest_inla <- function(
  fb,
  inla_fit,
  formula_slots,
  family_slots,
  reads,
  control_fixed,
  scalars = list(),
  completion = NULL,
  fixed_default_sd = numeric(0),
  route = "per_row"
) {
  inputs <- .fb_manifest_inputs(
    fb,
    scalars = scalars,
    completion = completion,
    fixed_default_sd = fixed_default_sd
  )
  readback <- .fb_manifest_inla_readback(inla_fit)
  rb <- readback[["rows"]]
  fixed_slots <- .fb_manifest_inla_fixed_slots(
    rb$tag[rb$section == "fixed"],
    control_fixed
  )
  emitted <- rbind(formula_slots, family_slots, fixed_slots)
  model_of <- stats::setNames(
    emitted$model,
    paste(emitted$section, emitted$tag)
  )
  group_model_of <- stats::setNames(
    emitted$group_model,
    paste(emitted$section, emitted$tag)
  )

  # A Wishart scale the package writes for a structure no stated prior
  # reaches is an input of its own.
  wishart <- emitted[
    !is.na(emitted$target_key) & startsWith(emitted$target_key, "wishart:"),
    ,
    drop = FALSE
  ]
  if (nrow(wishart)) {
    inputs <- rbind(inputs, .fb_manifest_input_frame(
      target = "wishart",
      target_key = wishart$tag,
      distribution = wishart$model,
      parameters = wishart$emitted_param,
      source = "structure_default"
    ))
  }

  # Every read-back row, matched to what was sent ---------------------------
  first_theta <- !duplicated(paste(rb$section, rb$tag))
  slot_rows <- list()
  matched <- logical(nrow(emitted))
  for (j in seq_len(nrow(rb))) {
    row <- rb[j, , drop = FALSE]
    candidates <- which(
      emitted$section == row$section & emitted$tag == row$tag
    )
    i <- if (length(candidates)) candidates[[1L]] else NA_integer_
    sent <- !is.na(i) && !isTRUE(emitted$pinned[[i]]) &&
      (!is.na(emitted$emitted_prior[[i]]) || !is.na(emitted$emitted_param[[i]]))
    # The hyperparameter a slot's prior lands on: the likelihood's residual
    # keyword, a fixed effect's own prior, or a latent term's first theta.
    carries_prior <- if (identical(row$section, "family")) {
      !is.na(i) &&
        identical(row$short_name, .fb_inla_hyper_keyword(emitted$model[[i]]))
    } else {
      identical(row$section, "fixed") || first_theta[[j]]
    }
    readback_text <- .fb_manifest_prior_text(row$prior, row$param, row$section)
    if (sent && carries_prior) {
      matched[[i]] <- TRUE
      slot_rows[[length(slot_rows) + 1L]] <- .fb_manifest_slot_out(
        section = row$section, tag = row$tag, theta = row$theta,
        target_key = emitted$target_key[[i]],
        input_id = .fb_manifest_input_for_key(
          row$section, emitted$target_key[[i]], inputs
        ),
        structure = emitted$structure[[i]],
        emitted = .fb_manifest_prior_text(
          emitted$emitted_prior[[i]], emitted$emitted_param[[i]], row$section
        ),
        readback = readback_text,
        agrees = .fb_manifest_same_prior(
          emitted$emitted_prior[[i]], emitted$emitted_param[[i]],
          row$prior, row$param
        ),
        source = "package",
        engine = "inla"
      )
      next
    }
    if (!is.na(i) && carries_prior) {
      matched[[i]] <- TRUE
    }
    model <- switch(
      row$section,
      fixed = row$tag,
      group = unname(group_model_of[paste("random", row$tag)]),
      unname(model_of[paste(row$section, row$tag)])
    )
    default <- if (is.na(model)) {
      NULL
    } else {
      .fb_manifest_inla_default(row$section, model, row$theta)
    }
    is_default <- !is.null(default) && .fb_manifest_same_prior(
      default$prior, .fb_manifest_param_string(default$param),
      row$prior, row$param
    )
    slot_rows[[length(slot_rows) + 1L]] <- .fb_manifest_slot_out(
      section = row$section, tag = row$tag, theta = row$theta,
      target_key = NA_character_, input_id = NA_character_,
      structure = if (is.na(i)) NA_character_ else emitted$structure[[i]],
      emitted = NA_character_, readback = readback_text,
      agrees = NA, source = if (is_default) "engine_default" else "unexpected",
      engine = "inla"
    )
  }

  # A prior the package sent that INLA does not report ----------------------
  lost <- which(
    !matched & (!is.na(emitted$emitted_prior) | !is.na(emitted$emitted_param)) &
      !emitted$pinned
  )
  for (i in lost) {
    slot_rows[[length(slot_rows) + 1L]] <- .fb_manifest_slot_out(
      section = emitted$section[[i]],
      tag = emitted$tag[[i]],
      theta = NA_character_,
      target_key = emitted$target_key[[i]],
      input_id = .fb_manifest_input_for_key(
        emitted$section[[i]], emitted$target_key[[i]], inputs
      ),
      structure = emitted$structure[[i]],
      emitted = .fb_manifest_prior_text(
        emitted$emitted_prior[[i]],
        emitted$emitted_param[[i]],
        emitted$section[[i]]
      ),
      readback = NA_character_, agrees = FALSE, source = "package",
      engine = "inla"
    )
  }
  slots <- .fb_manifest_bind_out(slot_rows)

  # Why an input reached no slot --------------------------------------------
  reasons <- c(
    "scalar:prior_vc_sd" = "prior_vc_sd_superseded_by_prior"
  )
  if ("sigma:sigma" %in% inputs$input_id &&
    !any(family_slots$target_key %in% "sigma")) {
    reasons[["sigma:sigma"]] <- "prior_family_has_no_residual_scale"
  }
  for (k in seq_len(nrow(reads))) {
    code <- switch(
      reads$structure[[k]],
      us_gxe = "us_gxe_prior_not_translatable_inla",
      simple_slope_cor = "cor_slope_prior_not_translatable_inla",
      NULL
    )
    if (!is.null(code)) {
      reasons[[paste0("sd:", reads$key[[k]])]] <- code
    }
  }

  .fb_manifest_assemble(
    engine = "inla",
    route = route,
    inputs = inputs,
    slots = slots,
    reasons = reasons,
    evidence = .fb_manifest_inla_evidence(inla_fit, route),
    readback_complete = readback[["complete"]] && !length(lost)
  )
}

.fb_manifest_slot_out <- function(
  section,
  tag,
  theta,
  target_key,
  input_id,
  structure,
  emitted,
  readback,
  agrees,
  source,
  engine
) {
  data.frame(
    slot_id = paste(section, tag, theta, sep = ":"),
    engine = engine,
    section = section,
    engine_slot = trimws(paste(section, tag, if (is.na(theta)) "" else theta)),
    target_key = target_key,
    input_id = input_id,
    structure = structure,
    emitted = emitted,
    readback = readback,
    readback_agrees = agrees,
    source = source,
    stringsAsFactors = FALSE
  )
}

.fb_manifest_bind_out <- function(rows) {
  template <- data.frame(
    slot_id = character(0),
    engine = character(0),
    section = character(0),
    engine_slot = character(0),
    target_key = character(0),
    input_id = character(0),
    structure = character(0),
    emitted = character(0),
    readback = character(0),
    readback_agrees = logical(0),
    source = character(0),
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, c(list(template), rows))
  rownames(out) <- NULL
  out
}

.fb_manifest_inla_evidence <- function(inla_fit, route) {
  if (route %in% c("aggregated_count", "streamed_count")) {
    return(.fb_manifest_evidence_unavailable("aggregated_evidence_unavailable"))
  }
  mlik <- inla_fit[["mlik"]]
  if (is.null(mlik) || !length(mlik) || !is.finite(mlik[1L, 1L])) {
    return(.fb_manifest_evidence_unavailable("evidence_not_computed"))
  }
  list(
    status = "available",
    value = unname(mlik[1L, 1L]),
    scale = "log marginal likelihood of the rows",
    reason_code = NA_character_
  )
}

.fb_manifest_evidence_unavailable <- function(reason_code) {
  list(
    status = "unavailable",
    value = NA_real_,
    scale = NA_character_,
    reason_code = reason_code
  )
}

#' Count what consumed each input, and put the manifest together
#'
#' @param reasons A named character vector: for an input that reached no
#'   slot, the reason code that explains it, keyed by input identifier.
#' @returns An `fb_manifest`.
#' @noRd
#' @keywords internal
.fb_manifest_assemble <- function(
  engine,
  route,
  inputs,
  slots,
  reasons,
  evidence,
  readback_complete
) {
  consumed <- slots[slots$source == "package" & !is.na(slots$input_id), ,
    drop = FALSE]
  n <- vapply(
    inputs$input_id,
    function(id) sum(consumed$input_id == id),
    integer(1L)
  )
  inputs$n_slots <- unname(n)
  inputs$cardinality <- ifelse(
    n == 0L,
    "not_consumed",
    ifelse(n == 1L, "one_to_one", "one_to_many")
  )
  inputs$expansion <- vapply(
    seq_len(nrow(inputs)),
    function(i) {
      if (n[[i]] < 2L) {
        return(NA_character_)
      }
      in_input <- consumed$input_id == inputs$input_id[[i]]
      structures <- unique(consumed$structure[in_input])
      sprintf("%d slots (%s)", n[[i]], paste(structures, collapse = ", "))
    },
    character(1L)
  )
  inputs$reason_code <- ifelse(
    n == 0L,
    unname(reasons[inputs$input_id]),
    NA_character_
  )

  structure(
    list(
      version = 1L,
      engine = engine,
      route = list(route = route),
      inputs = inputs,
      slots = slots,
      evidence = evidence,
      readback_complete = isTRUE(readback_complete)
    ),
    class = c("fb_manifest", "list")
  )
}


# --- brms ------------------------------------------------------------------ #

.fb_manifest_brms_tuple <- function(prior) {
  cols <- c("class", "coef", "group", "resp", "dpar", "nlpar")
  parts <- lapply(cols, function(col) {
    v <- as.character(prior[[col]] %||% rep("", NROW(prior)))
    v[is.na(v)] <- ""
    v
  })
  do.call(paste, c(parts, list(sep = "|")))
}

#' The manifest of a brms fit
#'
#' @param fb The `fb_terms` object the fit was built from.
#' @param emitted_prior The brms prior rows sent, after any merge.
#' @param brmsfit The `brmsfit`; an unsampled one (`empty = TRUE`) carries
#'   the same prior table.
#' @param tagged_prior The rows `.priors_to_brms()` built, carrying the input
#'   identifier of each row in the `fb_input_ids` attribute.
#' @param completion The default slope priors the brms emit divided by their
#'   covariates' standard deviations (`.fb_default_covariate_scaled()`), or
#'   `NULL`.
#' @returns An `fb_manifest`.
#' @noRd
#' @keywords internal
.fb_manifest_brms <- function(
  fb,
  emitted_prior,
  brmsfit,
  tagged_prior = emitted_prior,
  completion = NULL
) {
  emitted <- as.data.frame(emitted_prior %||% data.frame())
  tagged <- as.data.frame(tagged_prior %||% data.frame())
  tag_ids <- attr(tagged_prior, "fb_input_ids") %||%
    rep(NA_character_, nrow(tagged))
  tag_map <- stats::setNames(tag_ids, .fb_manifest_brms_tuple(tagged))

  emitted_tuple <- if (nrow(emitted)) {
    .fb_manifest_brms_tuple(emitted)
  } else {
    character(0)
  }
  input_ids <- unname(tag_map[emitted_tuple])

  # Scalar arguments the route used ------------------------------------------
  scalars <- list()
  for (i in seq_len(nrow(emitted))) {
    id <- input_ids[[i]]
    if (!is.na(id) && startsWith(id, "scalar:")) {
      scalars[[sub("^scalar:", "", id)]] <- as.character(emitted$prior[[i]])
    }
  }
  inputs <- .fb_manifest_inputs(
    fb,
    scalars = scalars[!duplicated(names(scalars))],
    completion = completion
  )

  # Rows with no specification behind them: prior rows passed to
  # fb_refit(extra_prior = ), and the log-sigma prior a sectioned residual
  # needs when nothing states one.
  extra_rows <- list()
  for (i in seq_len(nrow(emitted))) {
    id <- input_ids[[i]]
    if (!is.na(id) && id %in% inputs$input_id) {
      next
    }
    source <- if (is.na(id)) "extra_prior" else "structure_default"
    if (is.na(id)) {
      id <- paste0("extra_prior:", emitted_tuple[[i]])
      input_ids[[i]] <- id
    }
    density <- as.character(emitted$prior[[i]])
    extra_rows[[length(extra_rows) + 1L]] <- .fb_manifest_input_frame(
      target = sub(":.*$", "", id),
      target_key = sub("^[^:]*:", "", id),
      distribution = sub("\\(.*$", "", density),
      parameters = density,
      source = source
    )
  }
  if (length(extra_rows)) {
    extra <- do.call(rbind, extra_rows)
    inputs <- rbind(inputs, extra[!duplicated(extra$input_id), , drop = FALSE])
  }

  readback <- as.data.frame(brms::prior_summary(brmsfit))
  readback_tuple <- .fb_manifest_brms_tuple(readback)
  slot_id_of <- function(tuple) paste0("brms:", tuple)

  slot_rows <- list()
  found <- logical(nrow(emitted))
  for (i in seq_len(nrow(emitted))) {
    hit <- which(
      readback_tuple == emitted_tuple[[i]] & readback$source == "user"
    )
    found[[i]] <- length(hit) > 0L
    slot_rows[[length(slot_rows) + 1L]] <- data.frame(
      slot_id = slot_id_of(emitted_tuple[[i]]),
      engine = "brms",
      section = as.character(emitted$class[[i]]),
      engine_slot = emitted_tuple[[i]],
      target_key = emitted_tuple[[i]],
      input_id = input_ids[[i]],
      structure = as.character(emitted$class[[i]]),
      emitted = as.character(emitted$prior[[i]]),
      readback = if (found[[i]]) {
        as.character(readback$prior[[hit[[1L]]]])
      } else {
        NA_character_
      },
      readback_agrees = found[[i]] &&
        identical(trimws(as.character(emitted$prior[[i]])),
                  trimws(as.character(readback$prior[[hit[[1L]]]]))),
      source = "package",
      stringsAsFactors = FALSE
    )
  }
  for (j in seq_len(nrow(readback))) {
    if (
      readback$source[[j]] == "user" &&
        readback_tuple[[j]] %in% emitted_tuple
    ) {
      next
    }
    source <- switch(
      as.character(readback$source[[j]]),
      default = "engine_default",
      vectorized = "engine_expansion",
      "unexpected"
    )
    slot_rows[[length(slot_rows) + 1L]] <- data.frame(
      slot_id = slot_id_of(readback_tuple[[j]]),
      engine = "brms",
      section = as.character(readback$class[[j]]),
      engine_slot = readback_tuple[[j]],
      target_key = NA_character_,
      input_id = NA_character_,
      structure = as.character(readback$class[[j]]),
      emitted = NA_character_,
      readback = as.character(readback$prior[[j]]),
      readback_agrees = NA,
      source = source,
      stringsAsFactors = FALSE
    )
  }
  slots <- .fb_manifest_bind_out(slot_rows)

  reasons <- c(
    "scalar:prior_vc_sd" = "prior_vc_sd_superseded_by_prior"
  )
  if (
    "sigma:sigma" %in% inputs$input_id &&
      !.brms_family_has_sigma(fb[["family"]])
  ) {
    reasons[["sigma:sigma"]] <- "prior_family_has_no_residual_scale"
  }

  .fb_manifest_assemble(
    engine = "brms",
    route = "per_row",
    inputs = inputs,
    slots = slots,
    reasons = reasons,
    evidence = .fb_manifest_evidence_unavailable("evidence_not_computed"),
    readback_complete = all(found)
  )
}


# --- keys ------------------------------------------------------------------ #

#' Refuse a stated prior key that no engine slot reads
#'
#' @param keys The keys the caller's specifications are written under.
#' @param available The keys the emitter's slots read.
#' @param engine The engine name, for the message.
#' @returns `TRUE`, invisibly, when every key is read.
#' @noRd
#' @keywords internal
.fb_check_prior_keys_reachable <- function(keys, available, engine) {
  missing_keys <- setdiff(keys, available)
  if (!length(missing_keys)) {
    return(invisible(TRUE))
  }
  .fb_stop_prior_target_absent(
    target_label = paste0("sd(group = \"", missing_keys[[1L]], "\")"),
    kind = "variance component",
    available = available,
    engine = engine
  )
}
