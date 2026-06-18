# helper-genomic.R -- known-truth simulators for the genomics / MET
# parameter-recovery cells (G0c).
#
# Every recovery cell grounds the *config* as well as the data (the
# lasrosas Population=6 lesson): the simulator fixes the data-generating
# variance components / QTL effects so a fit's recovery is checked
# against truth the simulator authored, not against a self-consistent
# re-derivation. The relationship matrix is scaled to mean-diagonal 1
# (VanRaden convention) so the heritability target is unambiguous on the
# genotype-mean basis.

# Scale a relationship matrix to mean diagonal 1 so sigma_g^2 is the
# genetic variance on the genotype-mean basis and
# h^2 = sigma_g^2 / (sigma_g^2 + sigma_e^2) exactly.
scale_kinship_unit_diag <- function(K) {
  d_bar <- mean(diag(K))
  if (d_bar <= 0) {
    stop("kinship matrix has non-positive mean diagonal; cannot scale.")
  }
  K / d_bar
}

# A random VanRaden-style genomic relationship matrix from simulated
# biallelic markers, scaled to mean-diagonal 1 and ridged to strict PD.
sim_kinship <- function(n_geno = 25L, n_markers = 200L, maf = 0.25, seed = 1L) {
  set.seed(seed)
  Z <- matrix(stats::rbinom(n_geno * n_markers, 2L, maf), n_geno, n_markers)
  Zc <- scale(Z, center = TRUE, scale = FALSE)
  G <- tcrossprod(Zc) / n_markers
  G <- scale_kinship_unit_diag(G + diag(n_geno) * 1e-6)
  rownames(G) <- colnames(G) <- paste0("g", seq_len(n_geno))
  G
}

# Markers with founder relatedness: each genotype is the mid-parent of two
# random founders plus a little mutation, so genotypes share large genomic
# segments and the relationship matrix has real off-diagonal structure.
# Random independent markers (sim_kinship) make genotypes ~unrelated, where
# genomic prediction of a held-out line has no relatives to borrow from;
# this generator is the realistic genomic-selection setting.
sim_related_markers <- function(
  n_geno = 150L,
  n_markers = 600L,
  n_founders = 12L,
  maf = 0.25,
  mutation = 0.02,
  seed = 1L
) {
  set.seed(seed)
  founders <- matrix(
    stats::rbinom(n_founders * n_markers, 2L, maf), n_founders, n_markers
  )
  Z <- matrix(0L, n_geno, n_markers)
  for (i in seq_len(n_geno)) {
    par <- sample(n_founders, 2L)
    g <- round((founders[par[1L], ] + founders[par[2L], ]) / 2)
    flip <- stats::rbinom(n_markers, 1L, mutation)
    step <- flip * sample(c(-1L, 1L), n_markers, replace = TRUE)
    Z[i, ] <- pmin(2L, pmax(0L, g + step))
  }
  Z
}

# Genomic relationship matrix from a marker panel (VanRaden), scaled to
# unit diagonal and labelled.
kinship_from_markers <- function(Z) {
  Zc <- scale(Z, center = TRUE, scale = FALSE)
  G <- tcrossprod(Zc) / ncol(Z)
  G <- scale_kinship_unit_diag(G + diag(nrow(Z)) * 1e-6)
  rownames(G) <- colnames(G) <- paste0("g", seq_len(nrow(Z)))
  G
}

# Simulate a GBLUP phenotype with known variance components.
#
#   y_{ij} = mu + u_i + eps_{ij},
#   u ~ N(0, var_g * K),  eps ~ N(0, var_e).
#
# Returns the data frame plus the authored truth (var_g, var_e, h2, the
# realised breeding values u). The genetic effect is drawn as u = B a
# with B B' = K (spectral square root) and a ~ N(0, var_g I), so u has
# the exact target covariance.
sim_gblup_pheno <- function(
  K,
  var_g = 1,
  var_e = 1,
  n_rep = 4L,
  mu = 10,
  seed = 1L
) {
  set.seed(seed)
  n_geno <- nrow(K)
  e <- eigen((K + t(K)) / 2, symmetric = TRUE)
  vals <- pmax(e$values, 0)
  B <- sweep(e$vectors, 2L, sqrt(vals), `*`)
  u <- as.vector(B %*% stats::rnorm(n_geno, 0, sqrt(var_g)))
  names(u) <- rownames(K)
  geno <- factor(rep(rownames(K), times = n_rep), levels = rownames(K))
  y <- mu + u[as.character(geno)] + stats::rnorm(n_geno * n_rep, 0, sqrt(var_e))
  list(
    data = data.frame(geno = geno, y = as.numeric(y), stringsAsFactors = FALSE),
    K = K,
    u_true = u,
    var_g = var_g,
    var_e = var_e,
    h2 = var_g / (var_g + var_e),
    mu = mu
  )
}

# Simulate a GWAS phenotype with known QTL.
#
#   y_i = mu + sum_q Zstd[i, qtl_q] * effect_q + poly_i + eps_i,
#
# where the markers at `qtl_idx` carry the named additive effects, an
# optional polygenic background poly ~ N(0, var_poly * G) spreads small
# effects over all markers, and eps ~ N(0, var_e). One record per
# genotype (the GWAS unit). Returns the marker matrix, the phenotype,
# and the authored truth (which markers are QTL and their effects).
sim_gwas_pheno <- function(
  n_geno = 150L,
  n_markers = 300L,
  qtl_idx = c(40L, 120L, 220L),
  qtl_effect = c(1.5, -1.2, 0.9),
  var_e = 1,
  var_poly = 0,
  maf = 0.25,
  mu = 0,
  seed = 1L
) {
  set.seed(seed)
  if (length(qtl_idx) != length(qtl_effect)) {
    stop("qtl_idx and qtl_effect must be the same length.")
  }
  if (length(qtl_idx) && max(qtl_idx) > n_markers) {
    stop("qtl_idx (max ", max(qtl_idx), ") exceeds n_markers (", n_markers,
      "); pass an in-bounds qtl_idx.")
  }
  Z <- matrix(stats::rbinom(n_geno * n_markers, 2L, maf), n_geno, n_markers)
  colnames(Z) <- paste0("snp", seq_len(n_markers))
  Zstd <- scale(Z)
  Zstd[!is.finite(Zstd)] <- 0
  signal <- as.vector(Zstd[, qtl_idx, drop = FALSE] %*% qtl_effect)

  poly <- rep(0, n_geno)
  if (var_poly > 0) {
    Zc <- scale(Z, center = TRUE, scale = FALSE)
    G <- tcrossprod(Zc) / n_markers
    G <- scale_kinship_unit_diag(G + diag(n_geno) * 1e-6)
    e <- eigen(G, symmetric = TRUE)
    B <- sweep(e$vectors, 2L, sqrt(pmax(e$values, 0)), `*`)
    poly <- as.vector(B %*% stats::rnorm(n_geno, 0, sqrt(var_poly)))
  }

  y <- mu + signal + poly + stats::rnorm(n_geno, 0, sqrt(var_e))
  list(
    markers = Z,
    y = as.numeric(y),
    qtl_idx = qtl_idx,
    qtl_effect = qtl_effect,
    var_e = var_e,
    var_poly = var_poly,
    mu = mu
  )
}
