# pooled_prior -- partially pooled per-level standard deviations.
#
# `sd(group = "f:g", by = "f") ~ pooled(scale = <prior>)` and
# `sigma(by = "f") ~ pooled(scale = <prior>)` give the per-level SDs of
# diag(f):g and of dsum(~ units | f) one log-normal distribution with an
# estimated location and an estimated scale (the pooling SD). Fitted on brms
# only; the INLA gate refuses it (pooled_prior_requires_brms).
#
# The brms lowering is route B1 of the 2026-10-08 model C spike:
#   residual  `sigma ~ 1 + (1 | f)`, a normal(location, 2) prior on the
#             log-sigma intercept and the scale prior on its group SD;
#   diag      a flat prior on the per-level SDs and, through `stanvars`,
#             `target += lognormal_lpdf(sd_K | pool_mean_g, pool_sd_g)`.

# Distributions a pooling SD can take: every SD-scale density the brms
# translation carries.
.fb_pooled_scale_families <- function() {
  setdiff(.fb_prior_supported_families(), c("lkj", "pooled"))
}

# pooled(scale = <prior>), the scale prior parsed as a distribution of its
# own.
.parse_prior_pooled <- function(expr, envir) {
  scale <- as.list(expr[-1])[["scale"]]
  if (is.null(scale)) {
    stop(.fb_refusal_condition(
      reason_code = "prior_argument_missing",
      message = paste0(
        "`pooled()` needs a `scale` prior for the pooling SD, for example ",
        "pooled(scale = half_normal(scale = 0.5)). Got: ",
        .fb_prior_deparse(expr), "."
      )
    ))
  }
  inner <- .parse_prior_distribution(scale, envir = envir)
  if (!inner$family %in% .fb_pooled_scale_families()) {
    stop(.fb_refusal_condition(
      reason_code = "prior_distribution_unknown",
      message = paste0(
        "The pooling SD in ", .fb_prior_deparse(expr), " is given a `",
        inner$family, "()` prior. Supported: ",
        paste(.fb_pooled_scale_families(), collapse = ", "), "."
      )
    ))
  }
  list(family = "pooled", args = list(scale = inner))
}

# `by =` and pooled() come together: `by` names the factor whose levels are
# pooled, and only pooled() reads it.
.check_prior_pooled_pairing <- function(spec, formula) {
  has_by <- !is.null(spec$target$by)
  is_pooled <- identical(spec$spec$family, "pooled")
  if (is_pooled && !has_by) {
    stop(.fb_refusal_condition(
      reason_code = "prior_target_argument_missing",
      message = paste0(
        "`", .fb_prior_deparse(formula), "` pools a set of SDs without ",
        "naming the factor whose levels share the distribution. Write ",
        "sd(group = \"<f>:<g>\", by = \"<f>\") for diag(<f>):<g> or ",
        "sigma(by = \"<f>\") for dsum(~ units | <f>)."
      )
    ))
  }
  if (has_by && !is_pooled) {
    stop(.fb_refusal_condition(
      reason_code = "prior_target_unsupported",
      message = paste0(
        "`by =` in `", .fb_prior_deparse(formula), "` names the factor ",
        "whose levels pooled() pools over, and the distribution is not ",
        "pooled(). Write pooled(scale = ...) on the right, or drop `by =`."
      )
    ))
  }
  invisible(NULL)
}

.fb_pooled_specs <- function(prior) {
  if (!inherits(prior, "fb_prior")) {
    return(list())
  }
  Filter(function(s) identical(s$spec$family, "pooled"), prior$specs)
}

# One block per pooled() specification, resolved against the model's terms.
# A specification naming a term the model does not have is refused.
.fb_pooled_blocks <- function(fb) {
  specs <- .fb_pooled_specs(fb$priors)
  blocks <- list()
  for (s in specs) {
    by <- s$target$by
    label <- paste0(.format_prior_target(s$target), " ~ pooled()")
    if (identical(s$target$type, "sd")) {
      diag_terms <- Filter(
        function(t) identical(t$type %||% "", "at_simple") && is.null(t$level),
        fb$random_terms %||% list()
      )
      spelled <- vapply(
        diag_terms,
        function(t) paste0(t$outer, ":", t$inner),
        character(1L)
      )
      hit <- which(
        spelled == s$target$group &
          vapply(diag_terms, function(t) identical(t$outer, by), logical(1L))
      )
      if (!length(hit)) {
        .fb_stop_prior_target_absent(
          target_label = .format_prior_target(s$target),
          kind = "per-level variance term",
          available = paste0(
            "group = \"", spelled, "\", by = \"",
            vapply(diag_terms, function(t) t$outer, character(1L)), "\""
          )[nzchar(spelled)],
          engine = "brms"
        )
      }
      blocks[[length(blocks) + 1L]] <- list(
        kind = "sd",
        by = by,
        group = s$target$group,
        brms_group = diag_terms[[hit[[1L]]]]$inner,
        scale = s$spec$args$scale,
        input_id = .fb_manifest_spec_input_id(s),
        label = label
      )
      next
    }
    sections <- Filter(
      function(t) identical(t$type %||% "", "at_units") && is.null(t$level),
      fb$residual_terms %||% list()
    )
    vars <- vapply(sections, function(t) t$var, character(1L))
    if (!by %in% vars) {
      .fb_stop_prior_target_absent(
        target_label = .format_prior_target(s$target),
        kind = "sectioned residual",
        available = paste0("by = \"", vars, "\"")[nzchar(vars)],
        engine = "brms"
      )
    }
    plain <- Filter(
      function(x) {
        identical(x$target$type, "sigma") && !identical(x$spec$family, "pooled")
      },
      fb$priors$specs
    )
    if (length(plain)) {
      stop(.fb_refusal_condition(
        reason_code = "prior_argument_duplicated",
        message = paste0(
          "The residual SDs carry two priors: ",
          .format_prior_target(s$target), " ~ pooled() and sigma ~ ",
          .format_prior_distribution(plain[[1L]]$spec), ". Keep one."
        )
      ))
    }
    blocks[[length(blocks) + 1L]] <- list(
      kind = "sigma",
      by = by,
      group = "sigma",
      brms_group = by,
      scale = s$spec$args$scale,
      input_id = .fb_manifest_spec_input_id(s),
      label = label
    )
  }
  blocks
}

.fb_pooled_label <- function(blocks) {
  by <- unique(vapply(blocks, function(b) b$by, character(1L)))
  paste0("per-", tolower(paste(by, collapse = "/")), " SDs partially pooled")
}

# backend_decision()'s reason, with the label added when the model pools.
.fb_pooled_reason <- function(reason, fb) {
  specs <- .fb_pooled_specs(fb$priors)
  if (!length(specs)) {
    return(reason)
  }
  by <- unique(vapply(specs, function(s) s$target$by, character(1L)))
  paste0(
    reason, " Prior: per-", tolower(paste(by, collapse = "/")),
    " SDs partially pooled, fitted on brms only."
  )
}

# INLA gate check: a pooled prior is fitted on brms only.
.lgm_check_pooled_prior <- function(fb) {
  rule_id <- "pooled_prior_requires_brms"
  specs <- .fb_pooled_specs(fb$priors)
  if (!length(specs)) {
    return(.lgm_pass(rule_id))
  }
  spellings <- vapply(
    specs,
    function(s) {
      paste0(
        "`", .format_prior_target(s$target), " ~ ",
        .format_prior_distribution(s$spec), "`"
      )
    },
    character(1L)
  )
  .lgm_fail(
    rule_id,
    paste0(
      "pooled() is fitted on brms only: INLA cannot put a shared prior on ",
      "separate variance parameters. This model asks for ",
      paste(spellings, collapse = " and "),
      ". Pass backend = \"brms\", or leave backend = \"auto\", which ",
      "routes there."
    ),
    diagnostic = paste0("pooled() specifications: ", length(specs)),
    reason_code = rule_id
  )
}

# The location prior's centre on the log-SD scale: log(s / 5), with s the
# default prior's scale, which is log(sd(y)) on an identity link.
.fb_pooled_location <- function(fb, data) {
  s <- .default_half_normal_scale(data, fb$response, fb$family, fb$link)$scale
  log(s / .FB_DEFAULT_SD_MULTIPLIER)
}

# A Stan identifier from a brms group name.
.fb_stan_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  sub("_+$", "", x)
}

# A brms density string, "normal(0, 0.5)", as the Stan statement on `par`.
.fb_stan_lpdf <- function(density, par) {
  sub("^([a-z_]+)\\((.*)\\)$", paste0("\\1_lpdf(", par, " | \\2)"), density)
}

# brms's index K for the per-level SD vector `sd_K` of a group, read off the
# program brms writes when that group alone carries a marker prior.
.FB_POOLED_MARKER <- "normal(0, 987654.321)"

.fb_brms_sd_index <- function(formula, data, data2, family, group) {
  code <- brms::make_stancode(
    formula,
    data = data,
    data2 = data2,
    family = family,
    prior = brms::prior_string(.FB_POOLED_MARKER, class = "sd", group = group),
    silent = 2L
  )
  hit <- regmatches(
    code,
    regexpr("normal_lpdf\\(sd_[0-9]+ \\| 0, 987654\\.321\\)", code)
  )
  name <- sub("^normal_lpdf\\((sd_[0-9]+) .*$", "\\1", hit)
  if (length(hit) != 1L || !.fb_brms_stancode_declares(code, name)) {
    stop(.fb_refusal_condition(
      reason_code = "brms_cannot_represent_term",
      message = paste0(
        "flexyBayes could not find brms's per-level SD vector for group `",
        group, "` in the Stan program brms ",
        as.character(utils::packageVersion("brms")),
        " writes, so the pooled prior has no parameter to attach to."
      ),
      term_type = "pooled_prior"
    ))
  }
  as.integer(sub("^sd_", "", name))
}

# TRUE when `code` declares the per-level SD vector `name` ("sd_3").
.fb_brms_stancode_declares <- function(code, name) {
  if (length(name) != 1L || !grepl("^sd_[0-9]+$", name)) {
    return(FALSE)
  }
  k <- sub("^sd_", "", name)
  grepl(
    paste0("vector<lower=0>[M_", k, "] ", name, ";"),
    paste(code, collapse = "\n"),
    fixed = TRUE
  )
}

# The brms lowering of the pooled blocks. `specs` are the prior rows already
# built; the result carries the rows to add, the stanvars, and the Stan
# statements the read-back looks for.
.fb_brms_pooled_lowering <- function(blocks, specs, fb, formula, data, data2,
                                     family) {
  if (!length(blocks)) {
    return(NULL)
  }
  location <- .fb_pooled_location(fb, data)
  rows <- list()
  add_row <- function(string, class, group = NA_character_,
                      dpar = NA_character_, lb = NA_real_, ub = NA_real_,
                      input) {
    rows[[length(rows) + 1L]] <<- list(
      string = string, class = class, coef = NA_character_, group = group,
      dpar = dpar, lb = lb, ub = ub, input = input
    )
  }
  for (b in Filter(function(b) identical(b$kind, "sigma"), blocks)) {
    dens <- .brms_vc_density(b$scale)
    add_row("normal(pool_location, 2)", "Intercept", dpar = "sigma",
            input = b$input_id)
    add_row(dens$string, "sd", group = b$by, dpar = "sigma",
            lb = dens$lb, ub = dens$ub, input = b$input_id)
  }

  sd_blocks <- Filter(function(b) identical(b$kind, "sd"), blocks)
  par_code <- character(0)
  model_code <- character(0)
  statements <- list()
  if (length(sd_blocks)) {
    # brms's class-wide sd prior applies to every group without its own
    # row, and an empty group row does not override it. The class row is
    # made flat, and each other group of the mean keeps brms's default.
    avail <- as.data.frame(brms::default_prior(
      formula, data = data, data2 = data2, family = family
    ))
    for (col in c("coef", "group", "dpar")) {
      avail[[col]] <- as.character(avail[[col]] %||% "")
      avail[[col]][is.na(avail[[col]])] <- ""
    }
    in_mean <- avail$class == "sd" & avail$dpar == "" & avail$coef == ""
    class_default <- avail$prior[in_mean & avail$group == ""][[1L]]
    stated <- vapply(
      specs,
      function(s) {
        if (identical(s$class, "sd") && is.na(s$coef %||% NA_character_) &&
              is.na(s$dpar %||% NA_character_)) s$group else NA_character_
      },
      character(1L)
    )
    pooled_groups <- vapply(sd_blocks, function(b) b$brms_group, character(1L))
    fill <- setdiff(
      unique(avail$group[in_mean & avail$group != ""]),
      c(stated, pooled_groups)
    )
    add_row("", "sd", input = sd_blocks[[1L]]$input_id)
    for (g in fill) {
      add_row(class_default, "sd", group = g, lb = 0,
              input = paste0("structure:sd_", g))
    }
    for (b in sd_blocks) {
      k <- .fb_brms_sd_index(formula, data, data2, family, b$brms_group)
      nm <- .fb_stan_name(b$brms_group)
      mu <- paste0("pool_mean_", nm)
      tau <- paste0("pool_sd_", nm)
      dens <- .brms_vc_density(b$scale)
      upper <- if (is.finite(dens$ub)) paste0(", upper=", .fmt_num(dens$ub))
      par_code <- c(
        par_code,
        paste0("real ", mu, ";"),
        paste0("real<lower=", .fmt_num(dens$lb), upper, "> ", tau, ";")
      )
      lines <- c(
        hierarchy = sprintf(
          "target += lognormal_lpdf(sd_%d | %s, %s);", k, mu, tau
        ),
        location = sprintf("target += normal_lpdf(%s | pool_location, 2);", mu),
        scale = paste0("target += ", .fb_stan_lpdf(dens$string, tau), ";")
      )
      model_code <- c(model_code, lines)
      statements[[length(statements) + 1L]] <- list(
        block = b, sd_name = paste0("sd_", k), lines = lines,
        pooling_sd = tau, location = mu
      )
    }
  }

  stanvars <- brms::stanvar(x = location, name = "pool_location")
  if (length(par_code)) {
    stanvars <- stanvars +
      brms::stanvar(scode = paste(par_code, collapse = "\n  "),
                    block = "parameters") +
      brms::stanvar(scode = paste(unname(model_code), collapse = "\n  "),
                    block = "model")
  }
  list(
    specs = rows,
    stanvars = stanvars,
    statements = statements,
    blocks = blocks,
    location = location
  )
}

# Manifest slots for the pooled statements, read back from the program brms
# compiled. A statement missing from it, or on an SD vector it does not
# declare, does not agree.
.fb_manifest_brms_pooled_slots <- function(pooled, stancode) {
  if (is.null(pooled) || !length(pooled$statements)) {
    return(list())
  }
  code_lines <- trimws(strsplit(paste(stancode, collapse = "\n"), "\n")[[1L]])
  out <- list()
  for (st in pooled$statements) {
    declared <- .fb_brms_stancode_declares(stancode, st$sd_name)
    for (nm in names(st$lines)) {
      found <- st$lines[[nm]] %in% code_lines
      out[[length(out) + 1L]] <- data.frame(
        slot_id = paste0("stan:", st$block$brms_group, ":", nm),
        engine = "brms",
        section = "sd",
        engine_slot = if (identical(nm, "hierarchy")) st$sd_name else
          if (identical(nm, "location")) st$location else st$pooling_sd,
        target_key = st$block$group,
        input_id = st$block$input_id,
        structure = "pooled",
        emitted = st$lines[[nm]],
        readback = if (found) st$lines[[nm]] else NA_character_,
        readback_agrees = found && declared,
        source = "package",
        stringsAsFactors = FALSE
      )
    }
  }
  out
}

# What the fit records: the label, the blocks and each pooling SD's
# posterior median.
.fb_brms_pooled_record <- function(pooled, draws_mat) {
  if (is.null(pooled)) {
    return(NULL)
  }
  cn <- colnames(draws_mat)
  med <- numeric(0)
  for (b in pooled$blocks) {
    col <- if (identical(b$kind, "sd")) {
      paste0("pool_sd_", .fb_stan_name(b$brms_group))
    } else {
      paste0("sd_", b$by, "__sigma_Intercept")
    }
    if (col %in% cn) {
      key <- if (identical(b$kind, "sd")) b$group else "residual"
      med[[key]] <- stats::median(draws_mat[, col])
    }
  }
  list(
    label = .fb_pooled_label(pooled$blocks),
    targets = vapply(pooled$blocks, function(b) b$label, character(1L)),
    pooling_sd = med,
    location = pooled$location
  )
}

# The `prior` cell of summary()$varcomp for a component of a pooled block:
# the pooled() specification on each per-level SD, its scale prior on the
# pooling SD, NA for any other component.
.fb_varcomp_pooled_prior <- function(object, cmp) {
  fb <- object$extras$fb_terms
  if (is.null(object$extras$pooled) || is.null(fb)) {
    return(NA_character_)
  }
  for (s in .fb_pooled_specs(fb$priors)) {
    by <- s$target$by
    scale_prior <- .fb_prior_spec_string(list(spec = s$spec$args$scale))
    if (identical(s$target$type, "sigma")) {
      if (identical(cmp, paste0("sd_", by, "__sigma_Intercept"))) {
        return(scale_prior)
      }
      if (startsWith(cmp, "sigma_")) {
        return(.fb_prior_spec_string(s))
      }
      next
    }
    inner <- sub(paste0("^", by, ":"), "", s$target$group)
    if (startsWith(cmp, paste0("sd_", inner, "__", by))) {
      return(.fb_prior_spec_string(s))
    }
  }
  NA_character_
}

.populate_refusal_registry_pooled_prior <- function() {
  .register_refusal(
    reason_code = "pooled_prior_requires_brms",
    description = paste0(
      "A pooled() prior, one log-normal prior with an estimated mean and SD ",
      "shared by a set of per-level SDs, was sent to INLA. INLA cannot put ",
      "a shared prior on separate variance parameters, so flexyBayes fits ",
      "pooled() on brms only."
    ),
    message_template = paste0(
      "pooled() is fitted on brms only. This model asks for %s. Pass ",
      "backend = \"brms\", or leave backend = \"auto\", which routes there."
    ),
    registered_in_adr = paste0(
      "cairn 2026-10-09-model-c-hierarchical-per-trial-priors-join"
    ),
    plan_field = "rejected_routes",
    since_version = "0.17.0"
  )
}
