# Genomic triangulation. triangulate() compares two fits parameter by
# parameter; genomic quantities (heritability, variance components,
# breeding values, GWAS hits) need a genomics-aware comparison because the
# backends name them differently and the breeding values must be matched
# by genotype. These verbs harmonise via the same engine-agnostic
# extraction genomic_summary() uses, so a GBLUP fit on greta, INLA, or
# brms triangulates against another -- and, crucially, against an external
# field-standard *lens* (a REML answer from sommer, a GWAS hit list from
# GEMMA / rrBLUP) supplied as a plain list. flexyBayes core stays lean: it
# never depends on the field tools; the companion (flexyBayesOrchestra)
# builds the koine oracle that produces those lenses.
#
# As with triangulate(), this measures inter-lens *agreement*, not
# correspondence: two lenses fit to the same data share any fabricated
# upstream data fact (Independent Oracle Principle), so the result carries
# the same shared-upstream caveat.

# --- the genomic lens contract ------------------------------------ #

# Normalise either a flexybayes GBLUP / pedigree fit or a generic genomic
# lens list to the comparison shape:
#   list(h2, var_g, var_e = each c(mean, sd, lo, hi); gebv = named numeric;
#        label = character(1)).
# A flexybayes fit is extracted via .fb_genomic_draws() (posterior draws);
# a list lens may carry full draws (numeric vectors) or a point estimate
# with a standard error (the REML case), both reduced to the 4-number
# summary so a Bayesian posterior and a REML point estimate compare on the
# same footing.
.fb_genomic_lens <- function(x, term = NULL, label = NULL) {
  if (inherits(x, "flexybayes") || inherits(x, "flexybayes_inla")) {
    rt <- x$extras$parse_info$random %||% list()
    vmped <- Filter(function(t) (t$type %||% "") %in% c("vm", "ped"), rt)
    if (!length(vmped)) {
      stop(
        "triangulate_genomic(): a fit was supplied with no vm() / ped() ",
        "relationship term.",
        call. = FALSE
      )
    }
    var <- if (is.null(term)) vmped[[1L]]$var else term
    d <- .fb_genomic_draws(x, var)
    h2_draws <- d$var_g_draws / (d$var_g_draws + d$var_e_draws)
    h2_draws[!is.finite(h2_draws)] <- NA_real_
    gebv <- if (!is.null(d$gebv_draws)) {
      stats::setNames(colMeans(d$gebv_draws), d$labels)
    } else {
      NULL
    }
    return(list(
      h2 = .fb_lens_summ(h2_draws),
      var_g = .fb_lens_summ(d$var_g_draws),
      var_e = .fb_lens_summ(d$var_e_draws),
      gebv = gebv,
      label = label %||% .triangulate_source(x)
    ))
  }

  if (is.list(x)) {
    norm <- function(v) {
      if (is.null(v)) {
        return(c(mean = NA_real_, sd = NA_real_, lo = NA_real_, hi = NA_real_))
      }
      if (length(v) > 1L && is.null(names(v))) {
        return(.fb_lens_summ(v)) # draws
      }
      # point + optional se: list(estimate=, se=) or named c(estimate, se).
      est <- v[["estimate"]] %||% v[[1L]]
      se <- v[["se"]] %||% (if (length(v) > 1L) v[[2L]] else NA_real_)
      c(
        mean = as.numeric(est), sd = as.numeric(se),
        lo = as.numeric(est) - 1.96 * (se %||% 0),
        hi = as.numeric(est) + 1.96 * (se %||% 0)
      )
    }
    gebv <- x$gebv
    if (!is.null(gebv) && is.null(names(gebv))) {
      stop("a genomic lens `gebv` must be a *named* numeric vector ",
        "(names = genotype labels).", call. = FALSE)
    }
    return(list(
      h2 = norm(x$h2),
      var_g = norm(x$var_g),
      var_e = norm(x$var_e),
      gebv = gebv,
      label = label %||% x$label %||% "lens"
    ))
  }

  stop(
    "triangulate_genomic(): each argument must be a flexybayes GBLUP / ",
    "pedigree fit or a genomic-lens list(h2, var_g, var_e, gebv, label).",
    call. = FALSE
  )
}

.fb_lens_summ <- function(draws) {
  q <- stats::quantile(draws, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  c(
    mean = mean(draws, na.rm = TRUE), sd = stats::sd(draws, na.rm = TRUE),
    lo = q[[1L]], hi = q[[2L]]
  )
}

# --- triangulate_genomic ------------------------------------------ #

#' Triangulate genomic model outputs across engines or against a field lens
#'
#' Compare two genomic analyses -- heritability, variance components, and
#' genomic estimated breeding values -- on a common footing. Each argument
#' is either a flexyBayes GBLUP / pedigree fit (greta, INLA, or brms) or a
#' generic *genomic lens* (`list(h2, var_g, var_e, gebv, label)`), the form
#' a field-standard oracle such as sommer's REML supplies. Cross-engine use
#' checks that the Bayesian backends agree; the lens form lets the koine
#' fourth opinion (REML / established tools) cross-check the Bayesian
#' answer -- the orchestra's signature value in a field with decades of
#' established methods.
#'
#' Like [triangulate()], this measures inter-lens *agreement*, not
#' correspondence: two lenses fit to the same data share any fabricated
#' upstream data fact, so the result carries a shared-upstream caveat
#' unless `data_independence = TRUE`.
#'
#' @param a,b A flexyBayes GBLUP / pedigree fit, or a genomic-lens list. A
#'   lens entry for `h2` / `var_g` / `var_e` may be a numeric vector of
#'   draws or `list(estimate, se)` (a REML point); `gebv` is a named
#'   numeric vector keyed by genotype.
#' @param term Optional grouping-factor name selecting the relationship
#'   term when a fit has more than one.
#' @param data_independence `TRUE` if the two lenses used independently
#'   sourced data (suppresses the caveat); `NA` (default) or `FALSE`
#'   attaches it.
#'
#' @return A `triangulate_genomic_result`: `components` (a data frame
#'   comparing heritability and the variance components -- value and
#'   interval per lens, the difference, and an interval-overlap flag),
#'   `gebv` (the Pearson / Spearman correlation of the matched breeding
#'   values and the number in common), the lens labels, and the caveat.
#' @seealso [triangulate()], [genomic_summary()], [triangulate_gwas()].
#' @examples
#' # A Bayesian posterior (draws) against a field-standard REML point lens.
#' set.seed(1)
#' bv <- c(2, 1, 0, -1, -2)
#' bayes <- list(
#'   h2 = rnorm(500, 0.5, 0.05),
#'   gebv = stats::setNames(bv + rnorm(5, 0, 0.2), paste0("g", 1:5))
#' )
#' reml <- list(
#'   h2 = list(estimate = 0.48, se = 0.04),
#'   gebv = stats::setNames(bv + rnorm(5, 0, 0.2), paste0("g", 1:5))
#' )
#' triangulate_genomic(bayes, reml)
#' @export
triangulate_genomic <- function(a, b, term = NULL, data_independence = NA) {
  la <- .fb_genomic_lens(a, term)
  lb <- .fb_genomic_lens(b, term)

  comp_row <- function(name, sa, sb) {
    overlap <- !(sa[["hi"]] < sb[["lo"]] || sb[["hi"]] < sa[["lo"]])
    data.frame(
      quantity = name,
      value_a = sa[["mean"]], value_b = sb[["mean"]],
      sd_a = sa[["sd"]], sd_b = sb[["sd"]],
      difference = sa[["mean"]] - sb[["mean"]],
      intervals_overlap = isTRUE(overlap),
      row.names = NULL, stringsAsFactors = FALSE
    )
  }
  components <- rbind(
    comp_row("heritability", la$h2, lb$h2),
    comp_row("genetic_variance", la$var_g, lb$var_g),
    comp_row("residual_variance", la$var_e, lb$var_e)
  )

  gebv <- .fb_triangulate_gebv(la$gebv, lb$gebv)

  structure(
    list(
      components = components,
      gebv = gebv,
      label_a = la$label, label_b = lb$label,
      data_independence = data_independence,
      shared_upstream_caveat = .fb_genomic_caveat(data_independence)
    ),
    class = c("triangulate_genomic_result", "list")
  )
}

.fb_triangulate_gebv <- function(gebv_a, gebv_b) {
  if (is.null(gebv_a) || is.null(gebv_b)) {
    return(list(n_common = 0L, pearson = NA_real_, spearman = NA_real_))
  }
  common <- intersect(names(gebv_a), names(gebv_b))
  if (length(common) < 3L) {
    return(list(
      n_common = length(common), pearson = NA_real_, spearman = NA_real_
    ))
  }
  va <- gebv_a[common]
  vb <- gebv_b[common]
  list(
    n_common = length(common),
    pearson = stats::cor(va, vb),
    spearman = stats::cor(va, vb, method = "spearman")
  )
}

.fb_genomic_caveat <- function(data_independence) {
  if (isTRUE(data_independence)) {
    return(NA_character_)
  }
  paste0(
    "triangulate_genomic measures inter-lens agreement, not correspondence: ",
    if (identical(data_independence, FALSE)) {
      "both lenses used the SAME data, so agreement is common-mode. "
    } else {
      "data independence was not declared. "
    },
    "A fabricated upstream data fact (e.g. a mis-scaled relationship matrix) ",
    "is shared by both and would not be caught by their agreement."
  )
}

#' @exportS3Method print triangulate_genomic_result
print.triangulate_genomic_result <- function(x, ...) {
  cat("<triangulate_genomic_result>  ", x$label_a, " vs ", x$label_b, "\n",
    sep = "")
  for (i in seq_len(nrow(x$components))) {
    r <- x$components[i, ]
    if (is.na(r$value_a) || is.na(r$value_b)) {
      cat("  ", format(r$quantity, width = 18L), ": (not supplied)\n", sep = "")
      next
    }
    cat(
      "  ", format(r$quantity, width = 18L), ": ",
      format(round(r$value_a, 3L)), " vs ", format(round(r$value_b, 3L)),
      "  (diff ", format(round(r$difference, 3L)), ", ",
      if (r$intervals_overlap) "intervals overlap" else "INTERVALS DISJOINT",
      ")\n",
      sep = ""
    )
  }
  if (!is.null(x$gebv) && x$gebv$n_common >= 3L) {
    cat(
      "  GEBVs (", x$gebv$n_common, " common): r = ",
      format(round(x$gebv$pearson, 3L)), " (Pearson), ",
      format(round(x$gebv$spearman, 3L)), " (Spearman)\n",
      sep = ""
    )
  }
  if (!is.na(x$shared_upstream_caveat)) {
    cat("  ! ", x$shared_upstream_caveat, "\n", sep = "")
  }
  invisible(x)
}

# --- triangulate_gwas --------------------------------------------- #

#' Triangulate two genome-wide association scans
#'
#' Compare two GWAS results -- typically a flexyBayes [fb_gwas()] scan
#' against a field-standard scan (GEMMA / rrBLUP, supplied through the
#' koine oracle) -- by the agreement that matters for a scan: do the same
#' loci come up? Reports the Jaccard overlap of the genome-wide-significant
#' marker sets, the overlap among the top markers, the correlation of the
#' marker effects on the markers in common, and each scan's genomic-control
#' inflation factor.
#'
#' @param a,b An `fb_gwas` object, or a GWAS lens
#'   `list(results = <data frame with marker, p_value, effect>, lambda_gc)`.
#' @param alpha Genome-wide significance threshold on the Bonferroni-
#'   adjusted p-value for the hit-set comparison (default 0.05).
#' @param top_k Size of the top-marker overlap comparison (default 10).
#'
#' @return A `triangulate_gwas_result`: `jaccard` (significant-set
#'   overlap), `n_sig_a` / `n_sig_b` / `n_sig_common`, `top_k_overlap`,
#'   `effect_correlation` (on common markers), and `lambda_gc_a` /
#'   `lambda_gc_b`.
#' @seealso [fb_gwas()], [triangulate_genomic()].
#' @examples
#' mk <- function(sig) data.frame(
#'   marker = paste0("snp", 1:6),
#'   p_value = ifelse(seq_len(6) %in% sig, 1e-9, 0.5),
#'   p_bonferroni = ifelse(seq_len(6) %in% sig, 6e-9, 1),
#'   effect = c(2, 0, 0, -1.5, 0, 0)
#' )
#' a <- list(results = mk(c(1, 4)), lambda_gc = 1.01)
#' b <- list(results = mk(c(1, 4)), lambda_gc = 0.99)
#' triangulate_gwas(a, b)
#' @export
triangulate_gwas <- function(a, b, alpha = 0.05, top_k = 10L) {
  ra <- .fb_gwas_results(a)
  rb <- .fb_gwas_results(b)

  sig <- function(r) {
    padj <- r$p_bonferroni %||% pmin(1, r$p_value * nrow(r))
    r$marker[padj < alpha]
  }
  sa <- sig(ra)
  sb <- sig(rb)
  inter <- intersect(sa, sb)
  uni <- union(sa, sb)
  jaccard <- if (length(uni) == 0L) NA_real_ else length(inter) / length(uni)

  top <- function(r) r$marker[order(r$p_value)][seq_len(min(top_k, nrow(r)))]
  top_overlap <- length(intersect(top(ra), top(rb)))

  common <- intersect(ra$marker, rb$marker)
  eff_cor <- if (length(common) >= 3L) {
    ea <- ra$effect[match(common, ra$marker)]
    eb <- rb$effect[match(common, rb$marker)]
    suppressWarnings(stats::cor(ea, eb, use = "complete.obs"))
  } else {
    NA_real_
  }

  structure(
    list(
      jaccard = jaccard,
      n_sig_a = length(sa), n_sig_b = length(sb), n_sig_common = length(inter),
      top_k = top_k, top_k_overlap = top_overlap,
      effect_correlation = eff_cor,
      lambda_gc_a = .fb_gwas_lambda(a), lambda_gc_b = .fb_gwas_lambda(b)
    ),
    class = c("triangulate_gwas_result", "list")
  )
}

.fb_gwas_results <- function(x) {
  if (inherits(x, "fb_gwas")) {
    return(x$results)
  }
  if (is.list(x) && !is.null(x$results)) {
    return(as.data.frame(x$results))
  }
  stop("triangulate_gwas(): each argument must be an fb_gwas object or a ",
    "list(results, lambda_gc).", call. = FALSE)
}

.fb_gwas_lambda <- function(x) {
  if (inherits(x, "fb_gwas")) x$lambda_gc else (x$lambda_gc %||% NA_real_)
}

#' @exportS3Method print triangulate_gwas_result
print.triangulate_gwas_result <- function(x, ...) {
  cat("<triangulate_gwas_result>\n")
  cat(
    "  significant-set Jaccard: ", format(round(x$jaccard, 3L)),
    "  (", x$n_sig_common, " common of ", x$n_sig_a, " / ", x$n_sig_b, ")\n",
    sep = ""
  )
  cat("  top-", x$top_k, " overlap: ", x$top_k_overlap, "/", x$top_k, "\n",
    sep = "")
  cat("  effect correlation (common markers): ",
    format(round(x$effect_correlation, 3L)), "\n", sep = "")
  cat("  lambda_GC: ", format(round(x$lambda_gc_a, 3L)), " vs ",
    format(round(x$lambda_gc_b, 3L)), "\n", sep = "")
  invisible(x)
}
