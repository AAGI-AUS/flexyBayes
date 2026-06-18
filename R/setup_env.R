# Environment setup for flexyBayes
# Populates the evaluation environment with all R objects the generated code
# will reference (ID vectors, level counts, known matrices, etc.)
# Not exported.

# Set up the evaluation environment with data objects
#
# @param ev Environment to populate
# @param fixed_info Parsed fixed formula info
# @param random_terms List of random term descriptors
# @param rcov_terms List of rcov term descriptors
# @param data data.frame
# @param known_matrices Named list of matrices
# @param weights Optional numeric weight vector
.setup_env <- function(
  ev,
  fixed_info,
  random_terms,
  rcov_terms,
  data,
  known_matrices,
  weights
) {
  N <- nrow(data)
  ev$N_atg <- N
  ev$dat_atg <- data

  # Response
  resp <- fixed_info$response
  ev$y_atg <- data[[resp]]

  # Fixed effect id vectors
  for (term in fixed_info$terms) {
    if (term$type == "factor") {
      f <- factor(data[[term$var]])
      assign(paste0(term$var, "_id"), as.integer(f), envir = ev)
      assign(paste0("n_", term$var), nlevels(f), envir = ev)
    } else if (term$type == "continuous") {
      if (term$var %in% names(data)) {
        assign(term$var, as.numeric(data[[term$var]]), envir = ev)
      }
    } else if (term$type == "interaction") {
      for (v in term$vars) {
        if (v %in% names(data)) {
          col <- data[[v]]
          if (is.factor(col) || is.character(col)) {
            f <- factor(col)
            assign(paste0(v, "_id"), as.integer(f), envir = ev)
            assign(paste0("n_", v), nlevels(f), envir = ev)
          } else {
            assign(v, as.numeric(col), envir = ev)
          }
        }
      }
    } else if (term$type == "factor_interaction") {
      tag <- paste(term$vars, collapse = "_x_")
      combo <- do.call(
        paste,
        c(lapply(term$vars, function(v) data[[v]]), list(sep = "_"))
      )
      f <- factor(combo)
      assign(paste0(tag, "_id"), as.integer(f), envir = ev)
      assign(paste0("n_", tag), nlevels(f), envir = ev)
    } else if (term$type == "factor_numeric_interaction") {
      # v0.2.6 -- factor:continuous indexed
      # slopes. The codegen emits `slope_dev_<fac>_<con>[<fac>_id]`
      # so we need <fac>_id, n_<fac>, and the continuous data column.
      # For the Option D fallback (per-observation padded lookup) we
      # also pre-compute `is_ref_obs` (the reference-row mask) and
      # `shifted_idx` (the reference rows map to slot 1 but are
      # masked out before the multiply -- see codegen.R for the math).
      fac <- term$factor
      con <- term$continuous
      tag <- paste0(fac, "_", con)
      f <- factor(data[[fac]])
      f_id <- as.integer(f)
      assign(paste0(fac, "_id"), f_id, envir = ev)
      assign(paste0("n_", fac), nlevels(f), envir = ev)
      if (con %in% names(data)) {
        assign(con, as.numeric(data[[con]]), envir = ev)
      }
      # Option D scratch vectors -- harmless under Option C since
      # they live in the env but are never referenced. Computed
      # unconditionally so toggling `flexyBayes.force_option_d` at
      # fit time does not require re-running setup_env().
      assign(paste0(tag, "_is_ref_obs"), as.integer(f_id == 1L), envir = ev)
      assign(
        paste0(tag, "_shifted_idx"),
        as.integer(pmax(f_id - 1L, 1L)),
        envir = ev
      )
    }
  }

  # Random and rcov id vectors / level counts / matrices
  all_terms <- c(random_terms, rcov_terms)
  for (term in all_terms) {
    .setup_term_env(ev, term, data, known_matrices)
  }

  # Weights
  if (!is.null(weights)) ev$wt_atg <- weights
}

# v0.3.8 audit fix helpers. Resolve the grouping factor's
# level count + levels from data for known-matrix alignment validation.
# Return NULL when the variable is not present in `data` (e.g.
# data-side fixture binds *_id directly), in which case the validator
# downgrades to structural-only checks per its NULL contract.
.term_var_level_count <- function(vname, data) {
  if (!vname %in% names(data)) {
    return(NULL)
  }
  nlevels(factor(data[[vname]]))
}
.term_var_levels <- function(vname, data) {
  if (!vname %in% names(data)) {
    return(NULL)
  }
  levels(factor(data[[vname]]))
}

# Set up environment objects for a single random/rcov term
.setup_term_env <- function(ev, term, data, known_matrices) {
  mk_id <- function(vname) {
    if (vname %in% names(data)) {
      f <- factor(data[[vname]])
      assign(paste0(vname, "_id"), as.integer(f), envir = ev)
      assign(paste0("n_", vname), nlevels(f), envir = ev)
    } else if (paste0(vname, "_id") %in% names(data)) {
      id_col <- data[[paste0(vname, "_id")]]
      assign(paste0(vname, "_id"), id_col, envir = ev)
      assign(paste0("n_", vname), max(id_col), envir = ev)
    } else {
      warning("Variable '", vname, "' not found in data.")
    }
  }
  mk_mat <- function(mat_name) {
    if (is.null(mat_name) || is.na(mat_name)) {
      return()
    }
    if (!mat_name %in% names(known_matrices)) {
      stop(
        "Matrix '",
        mat_name,
        "' not found in known_matrices. ",
        "Provide it via known_matrices = list(",
        mat_name,
        " = <your matrix>)."
      )
    }
    assign(mat_name, known_matrices[[mat_name]], envir = ev)
  }

  switch(
    term$type,
    "simple" = ,
    "ide" = ,
    "id" = mk_id(term$var),
    # Uncorrelated random slope: bind the grouping-factor
    # ID + level count via the standard mk_id() helper, AND bind the
    # slope variable as a numeric column the codegen can lift via
    # as_data().
    "simple_slope_uncor" = {
      mk_id(term$var)
      sv <- term$slope_var
      if (!is.null(sv) && nzchar(sv)) {
        if (sv %in% names(data)) {
          assign(sv, as.numeric(data[[sv]]), envir = ev)
        } else {
          stop(
            "Slope variable '",
            sv,
            "' not found in data ",
            "(random-slope term '(",
            sv,
            " || ",
            term$var,
            ")').",
            call. = FALSE
          )
        }
      }
    },
    "vm" = {
      # v0.3.8 audit fix: invoke mk_id() before the route
      # check so the validators can enforce dim + level alignment of
      # the user-supplied known matrix against the grouping factor.
      mk_id(term$var)
      .stage5a_route_check(
        term,
        known_matrices,
        expected_n = .term_var_level_count(term$var, data),
        fit_levels = .term_var_levels(term$var, data)
      )
      cov <- term$cov_representation
      mk_mat(if (!is.null(cov)) cov$data else term$mat)
    },
    "ped" = {
      mk_id(term$var)
      .stage5a_route_check(
        term,
        known_matrices,
        expected_n = .term_var_level_count(term$var, data),
        fit_levels = .term_var_levels(term$var, data)
      )
      cov <- term$cov_representation
      mk_mat(if (!is.null(cov)) cov$data else term$mat)
    },
    "fa" = mk_id(term$var),
    "fa_gxe" = {
      mk_id(term$outer)
      mk_id(term$inner)
    },
    "us_gxe" = {
      mk_id(term$outer)
      mk_id(term$inner)
    },
    "at_simple" = {
      mk_id(term$outer)
      mk_id(term$inner)
    },
    "at_units" = mk_id(term$var),
    "nested" = {
      tag <- paste0(term$inner, "_in_", term$outer)
      nid <- as.integer(factor(paste(
        data[[term$outer]],
        data[[term$inner]],
        sep = "_"
      )))
      ev[[paste0("nested_id_", tag)]] <- nid
      ev[[paste0("n_", tag)]] <- length(unique(nid))
    },
    "combo" = {
      tag <- paste(rev(term$vars), collapse = "_in_")
      nid <- as.integer(factor(do.call(
        paste,
        c(lapply(term$vars, function(v) data[[v]]), list(sep = "_"))
      )))
      ev[[paste0("combo_id_", tag)]] <- nid
      ev[[paste0("n_", tag)]] <- length(unique(nid))
    },
    "ar1_spatial" = {
      if (!is.null(term$row_var) && term$row_var %in% names(data)) {
        row_vals <- data[[term$row_var]]
        assign(
          paste0(term$row_var, "_id"),
          as.integer(factor(row_vals)),
          envir = ev
        )
        assign(paste0("n_", term$row_var), length(unique(row_vals)), envir = ev)
      }
      if (!is.null(term$col_var) && term$col_var %in% names(data)) {
        col_vals <- data[[term$col_var]]
        assign(
          paste0(term$col_var, "_id"),
          as.integer(factor(col_vals)),
          envir = ev
        )
        assign(paste0("n_", term$col_var), length(unique(col_vals)), envir = ev)
      }
    }
  )
}
