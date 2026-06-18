# Genomic output contract (G0b) --- the standardised summary every
# genomic / MET fit carries on `fit$extras$genomic`.
#
# `.fb_genomic_summary()` turns posterior draws into the four quantities
# a breeder reads off a genomic-selection fit:
#
#   * heritability   h^2 = sigma_g^2 / (sigma_g^2 + sigma_e^2), the
#                    first quantity of interest (narrow-sense, on the
#                    genotype-mean basis; kinship scaling is the
#                    caller's declared convention --- see the MET /
#                    genomics vignette pitfall).
#   * GEBVs          genomic estimated breeding values u_i with full
#                    posterior uncertainty (mean + credible interval),
#                    the ranking quantity selection acts on.
#   * reliability    r_i^2 = 1 - PEV_i / sigma_g^2, the per-genotype
#                    prediction reliability (PEV_i = posterior variance
#                    of the BLUP). Clamped to [0, 1].
#   * marker effects (optional) per-marker posterior summary + the
#                    posterior probability the effect is retained, for
#                    whole-genome marker regression (G4).
#
# It is engine-agnostic: it consumes draws (numeric vectors / matrices),
# never a backend object, so the greta / INLA / brms emit paths all feed
# the same constructor and the result is triangulatable. The fit-level
# accessor that extracts the draws from a fitted model and populates
# `fit$extras$genomic` lands with the GBLUP emit route (G1); this file
# ships the pure constructor + schema + display the recovery cells and
# emit paths build on.

# --- constructor -------------------------------------------------- #

# `.fb_genomic_summary()` --- build the genomic-extras object from
# posterior draws.
#
# Arguments
#   var_g_draws  Numeric vector of additive genetic variance (sigma_g^2)
#                draws. Required.
#   var_e_draws  Numeric vector of residual variance (sigma_e^2) draws,
#                same length. Required.
#   gebv_draws   Optional n_draws x n_genotypes matrix of breeding-value
#                (u) draws. When supplied, GEBVs + reliability are
#                computed.
#   marker_draws Optional n_draws x n_markers matrix of marker-effect
#                draws (G4). When supplied, marker effects + posterior
#                retention probabilities are computed.
#   labels       Optional character labels for the genotype columns of
#                `gebv_draws`. Defaults to the column names or g1..gK.
#   marker_labels
#                Optional character labels for the marker columns.
#   inclusion_eps
#                Threshold below which a marker effect counts as
#                effectively zero when estimating the posterior
#                retention probability (G4 spike-and-slab style).
#                Default 0 (any non-zero draw counts as retained), which
#                degrades to NA for continuous-shrinkage draws that are
#                never exactly zero --- the caller passes a small eps for
#                those.
#
# Returns an `fb_genomic_summary` object: a classed list with
#   $heritability, $genetic_variance, $residual_variance
#       named numeric c(mean, sd, q2.5, q97.5).
#   $gebv          data.frame(genotype, mean, sd, q2.5, q97.5,
#                  reliability) or NULL.
#   $marker_effects data.frame(marker, mean, sd, q2.5, q97.5,
#                  prob_retained) or NULL.
#   $n_draws, $n_genotypes, $n_markers metadata.
.fb_genomic_summary <- function(
  var_g_draws,
  var_e_draws,
  gebv_draws = NULL,
  marker_draws = NULL,
  labels = NULL,
  marker_labels = NULL,
  inclusion_eps = 0
) {
  .fb_genomic_check_draws(var_g_draws, "var_g_draws")
  .fb_genomic_check_draws(var_e_draws, "var_e_draws")
  if (length(var_g_draws) != length(var_e_draws)) {
    stop(
      "`var_g_draws` and `var_e_draws` must be the same length (one ",
      "value per posterior draw); got ", length(var_g_draws), " and ",
      length(var_e_draws), ".",
      call. = FALSE
    )
  }
  if (any(var_g_draws < 0) || any(var_e_draws < 0)) {
    stop(
      "variance draws must be non-negative; got a negative `var_g_draws` ",
      "or `var_e_draws` entry. Pass variances (sigma^2), not SDs that ",
      "were squared incorrectly, and check the draw extraction.",
      call. = FALSE
    )
  }
  n_draws <- length(var_g_draws)

  h2_draws <- var_g_draws / (var_g_draws + var_e_draws)
  # A draw with sigma_g^2 = sigma_e^2 = 0 is degenerate (0/0); treat its
  # heritability as NA rather than NaN so the summary is honest.
  h2_draws[!is.finite(h2_draws)] <- NA_real_

  out <- list(
    heritability = .fb_genomic_summ(h2_draws),
    genetic_variance = .fb_genomic_summ(var_g_draws),
    residual_variance = .fb_genomic_summ(var_e_draws),
    gebv = NULL,
    marker_effects = NULL,
    n_draws = n_draws,
    n_genotypes = 0L,
    n_markers = 0L
  )

  if (!is.null(gebv_draws)) {
    gebv_draws <- .fb_genomic_as_draw_matrix(gebv_draws, n_draws, "gebv_draws")
    n_geno <- ncol(gebv_draws)
    labels <- .fb_genomic_labels(labels, gebv_draws, n_geno, "g")
    var_g_hat <- mean(var_g_draws)
    pev <- apply(gebv_draws, 2L, stats::var)
    reliability <- if (var_g_hat > 0) {
      pmax(0, pmin(1, 1 - pev / var_g_hat))
    } else {
      rep(NA_real_, n_geno)
    }
    out$gebv <- data.frame(
      genotype = labels,
      mean = colMeans(gebv_draws),
      sd = apply(gebv_draws, 2L, stats::sd),
      q2.5 = apply(gebv_draws, 2L, stats::quantile, 0.025, names = FALSE),
      q97.5 = apply(gebv_draws, 2L, stats::quantile, 0.975, names = FALSE),
      reliability = reliability,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    out$n_genotypes <- n_geno
  }

  if (!is.null(marker_draws)) {
    marker_draws <- .fb_genomic_as_draw_matrix(
      marker_draws, n_draws, "marker_draws"
    )
    n_markers <- ncol(marker_draws)
    marker_labels <- .fb_genomic_labels(
      marker_labels, marker_draws, n_markers, "m"
    )
    prob_retained <- if (inclusion_eps > 0) {
      colMeans(abs(marker_draws) > inclusion_eps)
    } else {
      rep(NA_real_, n_markers)
    }
    out$marker_effects <- data.frame(
      marker = marker_labels,
      mean = colMeans(marker_draws),
      sd = apply(marker_draws, 2L, stats::sd),
      q2.5 = apply(marker_draws, 2L, stats::quantile, 0.025, names = FALSE),
      q97.5 = apply(marker_draws, 2L, stats::quantile, 0.975, names = FALSE),
      prob_retained = prob_retained,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    out$n_markers <- n_markers
  }

  structure(out, class = c("fb_genomic_summary", "list"))
}

# --- display ------------------------------------------------------ #

#' @exportS3Method print fb_genomic_summary
print.fb_genomic_summary <- function(x, ...) {
  cat("<fb_genomic_summary>  (", x$n_draws, " draws)\n", sep = "")
  fmt <- function(s) {
    paste0(
      format(round(s[["mean"]], 4L)), " [",
      format(round(s[["q2.5"]], 4L)), ", ",
      format(round(s[["q97.5"]], 4L)), "]"
    )
  }
  cat("  heritability h^2 : ", fmt(x$heritability), "\n", sep = "")
  cat("  genetic variance : ", fmt(x$genetic_variance), "\n", sep = "")
  cat("  residual variance: ", fmt(x$residual_variance), "\n", sep = "")
  if (!is.null(x$gebv)) {
    cat(
      "  GEBVs            : ", x$n_genotypes, " genotypes (mean reliability ",
      format(round(mean(x$gebv$reliability, na.rm = TRUE), 3L)), ")\n",
      sep = ""
    )
  }
  if (!is.null(x$marker_effects)) {
    cat("  marker effects   : ", x$n_markers, " markers\n", sep = "")
  }
  invisible(x)
}

# --- fit-level accessor (G1) -------------------------------------- #

#' Genomic summary of a fitted relationship model
#'
#' Extract the breeder-facing genomic quantities -- narrow-sense
#' heritability \eqn{h^2}, genomic estimated breeding values (GEBVs) with
#' posterior reliability, and the genetic / residual variances -- from a
#' fitted `vm()` (genomic / GBLUP) or `ped()` (pedigree) model. The
#' quantities are read from the posterior draws engine-agnostically: a
#' greta, INLA, or brms GBLUP fit returns the same summary object, so a
#' multi-backend genomic analysis is directly triangulatable.
#'
#' The heritability is computed per draw as
#' \eqn{h^2 = \sigma_g^2 / (\sigma_g^2 + \sigma_e^2)} on the
#' genotype-mean basis; the kinship scaling convention is the analyst's
#' (state it when reporting). Reliability is
#' \eqn{1 - \mathrm{PEV}_i / \sigma_g^2} from the posterior variance of
#' each breeding value. GEBVs are available on the brms and INLA backends
#' natively and on the greta backend (the breeding-value vector is
#' monitored).
#'
#' @param fit A fitted `flexybayes` object carrying at least one `vm()`
#'   or `ped()` relationship term.
#' @param term Optional grouping-factor name selecting which relationship
#'   term to summarise when the model has more than one. Defaults to the
#'   first.
#'
#' @return An `fb_genomic_summary` object: `heritability`,
#'   `genetic_variance`, `residual_variance` (each a posterior summary),
#'   `gebv` (data frame of breeding values with reliability), and
#'   metadata.
#'
#' @seealso [fb_structured_cov()] for factor-analytic MET covariance,
#'   [triangulate()] for cross-engine agreement.
#' @examples
#' \dontrun{
#' fit <- flexybayes(
#'   yield ~ 1, random = ~ vm(geno, Gmat), data = met,
#'   known_matrices = list(Gmat = G), backend = "greta"
#' )
#' gs <- genomic_summary(fit)
#' gs$heritability
#' head(gs$gebv)
#' }
#' @export
genomic_summary <- function(fit, term = NULL) {
  if (!inherits(fit, "flexybayes") && !inherits(fit, "flexybayes_inla")) {
    stop("`fit` must be a flexybayes object.", call. = FALSE)
  }
  rt <- fit$extras$parse_info$random %||% list()
  vmped <- Filter(function(t) (t$type %||% "") %in% c("vm", "ped"), rt)
  if (!length(vmped)) {
    stop(
      "genomic_summary(): this fit carries no vm() / ped() relationship ",
      "term, so there is no genomic quantity to summarise. Fit a GBLUP / ",
      "pedigree model, e.g. random = ~ vm(geno, Gmat).",
      call. = FALSE
    )
  }
  chosen <- if (is.null(term)) {
    if (length(vmped) > 1L) {
      message(
        "flexyBayes: genomic_summary() found ", length(vmped),
        " relationship terms; summarising the first ('", vmped[[1L]]$var,
        "'). Pass term = to choose another."
      )
    }
    vmped[[1L]]
  } else {
    hit <- Filter(function(t) identical(t$var, term), vmped)
    if (!length(hit)) {
      stop(
        "genomic_summary(): '", term, "' is not a vm() / ped() term in ",
        "this fit. Available: ",
        paste(vapply(vmped, function(t) t$var, character(1)), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    hit[[1L]]
  }

  d <- .fb_genomic_draws(fit, chosen$var)
  .fb_genomic_summary(
    d$var_g_draws, d$var_e_draws,
    gebv_draws = d$gebv_draws, labels = d$labels
  )
}

# Backend-aware extraction of the genetic / residual variance draws and
# the breeding-value draw matrix from a fitted relationship model. Each
# backend names these quantities differently; this is the one place that
# knows the mapping (greta: sigma_<var> / sigma_e_atg / u_<var>[i,1];
# brms: sd_<var>__Intercept / sigma / r_<var>[lev,Intercept]; INLA:
# Precision for <var>_id / Precision for the Gaussian observations /
# <var>_id:i). Variances are returned (SDs squared, precisions inverted).
.fb_genomic_draws <- function(fit, var) {
  dr <- fb_as_draws_simple(fit)
  nm <- names(dr)
  # Canonical genotype labels = the relationship-matrix dimnames, which the
  # known-matrix alignment contract guarantees equal levels(<group>) in
  # order. Using them makes GEBV labels identical across greta / INLA /
  # brms, so triangulate_genomic() can match breeding values by genotype.
  lv <- .fb_genomic_levels(fit, var)

  collect <- function(cols) {
    if (!length(cols)) {
      return(NULL)
    }
    as.matrix(do.call(cbind, dr[cols]))
  }
  label_by_level <- function(gcols, fallback) {
    if (!is.null(lv) && length(lv) == length(gcols)) lv else fallback
  }

  if (inherits(fit, "flexybayes_inla")) {
    sg <- dr[[paste0("Precision for ", var, "_id")]]
    se <- dr[["Precision for the Gaussian observations"]]
    var_g <- if (!is.null(sg)) 1 / sg else NULL
    var_e <- if (!is.null(se)) 1 / se else NULL
    gcols <- grep(paste0("^", var, "_id:[0-9]+$"), nm, value = TRUE)
    gcols <- gcols[order(as.integer(sub(".*:", "", gcols)))]
    labels <- label_by_level(gcols, sub(paste0("^", var, "_id:"), "", gcols))
  } else if (inherits(fit, "flexybayes_brms")) {
    sg <- dr[[paste0("sd_", var, "__Intercept")]]
    se <- dr[["sigma"]]
    var_g <- if (!is.null(sg)) sg^2 else NULL
    var_e <- if (!is.null(se)) se^2 else NULL
    gcols <- grep(paste0("^r_", var, "\\["), nm, value = TRUE)
    # brms already names breeding values by the factor level.
    labels <- sub(paste0("^r_", var, "\\[(.*),Intercept\\]$"), "\\1", gcols)
  } else {
    sg <- dr[[paste0("sigma_", var)]]
    se <- dr[["sigma_e_atg"]]
    var_g <- if (!is.null(sg)) sg^2 else NULL
    var_e <- if (!is.null(se)) se^2 else NULL
    gcols <- grep(paste0("^u_", var, "\\["), nm, value = TRUE)
    gcols <- gcols[order(as.integer(sub("^.*\\[([0-9]+).*", "\\1", gcols)))]
    labels <- label_by_level(gcols, paste0(var, seq_along(gcols)))
  }

  if (is.null(var_g) || is.null(var_e)) {
    stop(
      "genomic_summary(): could not locate the genetic and residual ",
      "variance draws for term '", var, "' on the ",
      class(fit)[[1L]], " backend. The fit may predate the genomic ",
      "output contract, or the term name does not match the draws.",
      call. = FALSE
    )
  }

  list(
    var_g_draws = as.numeric(var_g),
    var_e_draws = as.numeric(var_e),
    gebv_draws = collect(gcols),
    labels = if (length(labels)) labels else NULL
  )
}

# --- internal helpers --------------------------------------------- #

# Canonical genotype labels for a vm() / ped() term: the dimnames of the
# relationship matrix it carries. The known-matrix alignment contract
# (.check_known_matrix_dimnames) guarantees these equal levels(<group>) in
# fit order, so they label breeding values identically on every backend.
# Returns NULL when the matrix has no dimnames (the labels fall back to the
# backend-native form).
.fb_genomic_levels <- function(fit, var) {
  rt <- fit$extras$parse_info$random %||% list()
  term <- NULL
  for (t in rt) {
    if (identical(t$var, var) && (t$type %||% "") %in% c("vm", "ped")) {
      term <- t
      break
    }
  }
  if (is.null(term)) {
    return(NULL)
  }
  sym <- (term$cov_representation$data %||% term$mat)
  km <- fit$extras$call_info$known_matrices
  if (is.null(km) || is.null(sym) || is.na(sym)) {
    return(NULL)
  }
  m <- km[[sym]]
  if (is.null(m)) {
    return(NULL)
  }
  rn <- rownames(m)
  if (is.null(rn)) NULL else as.character(rn)
}

# Five-number posterior summary used across the genomic quantities.
.fb_genomic_summ <- function(draws) {
  q <- stats::quantile(draws, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  c(
    mean = mean(draws, na.rm = TRUE),
    sd = stats::sd(draws, na.rm = TRUE),
    q2.5 = q[[1L]],
    q97.5 = q[[2L]]
  )
}

# Validate a draw vector: non-empty finite numeric.
.fb_genomic_check_draws <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L) {
    stop("`", name, "` must be a non-empty numeric vector of draws.",
      call. = FALSE
    )
  }
  if (any(!is.finite(x))) {
    stop("`", name, "` contains non-finite draws (NA / NaN / Inf).",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Coerce a draw matrix and check its row count matches n_draws.
.fb_genomic_as_draw_matrix <- function(M, n_draws, name) {
  M <- as.matrix(M)
  if (!is.numeric(M)) {
    stop("`", name, "` must be a numeric draws matrix.", call. = FALSE)
  }
  if (nrow(M) != n_draws) {
    stop(
      "`", name, "` must have one row per posterior draw (", n_draws,
      "); got ", nrow(M), " rows.",
      call. = FALSE
    )
  }
  if (any(!is.finite(M))) {
    stop("`", name, "` contains non-finite draws.", call. = FALSE)
  }
  M
}

# Resolve column labels from an explicit vector, the matrix colnames, or
# a generated prefix.
.fb_genomic_labels <- function(labels, M, n, prefix) {
  if (!is.null(labels)) {
    if (length(labels) != n) {
      stop(
        "label vector length (", length(labels), ") does not match the ",
        "number of columns (", n, ").",
        call. = FALSE
      )
    }
    return(as.character(labels))
  }
  if (!is.null(colnames(M))) {
    return(colnames(M))
  }
  paste0(prefix, seq_len(n))
}
