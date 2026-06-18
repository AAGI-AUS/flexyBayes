# Generate the met_example dataset for flexyBayes
# Run this script to recreate data/met_example.rda

set.seed(42)

n_geno <- 10
n_env <- 4
n_rep <- 3
N <- n_geno * n_env * n_rep # 120

# Base experimental design
dat <- expand.grid(
  geno = factor(paste0("G", seq_len(n_geno))),
  env = factor(paste0("E", seq_len(n_env))),
  rep = factor(paste0("R", seq_len(n_rep)))
)

# Spatial layout
dat$block <- factor(paste0("B", rep(seq_len(6), length.out = N)))
dat$row <- factor(rep(seq_len(10), length.out = N))
dat$col <- factor(rep(seq_len(12), length.out = N))

# True effects
mu <- 50
env_eff <- c(-3, -1, 1, 3) # environment effects
gen_eff <- rnorm(n_geno, 0, 2) # genotype effects (sigma_g = 2)
gxe_eff <- matrix(rnorm(n_geno * n_env, 0, sqrt(2)), n_geno, n_env)

# Generate yield
dat$yield <- mu +
  env_eff[as.integer(dat$env)] +
  gen_eff[as.integer(dat$geno)] +
  gxe_eff[cbind(as.integer(dat$geno), as.integer(dat$env))] +
  rnorm(N, 0, 1)

# Continuous covariate
dat$x_cov <- rnorm(N)

# Binary response
p <- plogis(-1 + 0.5 * scale(dat$yield))
dat$bin_y <- rbinom(N, 1, p)

# Count response
lambda <- exp(1 + 0.02 * dat$yield)
dat$count_y <- rpois(N, lambda)

# Genomic marker matrix (50 SNPs)
M_geno <- matrix(
  sample(c(0, 1, 2), n_geno * 50, replace = TRUE, prob = c(0.25, 0.5, 0.25)),
  nrow = n_geno,
  ncol = 50
)
rownames(M_geno) <- levels(dat$geno)
colnames(M_geno) <- paste0("SNP", seq_len(50))

# Genomic relationship matrix (VanRaden method 1)
p_freq <- colMeans(M_geno) / 2
W <- scale(M_geno, center = 2 * p_freq, scale = FALSE)
denom <- 2 * sum(p_freq * (1 - p_freq))
G_mat <- (W %*% t(W)) / denom
# Ensure positive definite
G_mat <- G_mat + diag(0.01, n_geno)
rownames(G_mat) <- colnames(G_mat) <- levels(dat$geno)

# Simple pedigree matrix (half-sib structure)
A_mat <- diag(n_geno)
for (i in 1:(n_geno - 1)) {
  for (j in (i + 1):n_geno) {
    A_mat[i, j] <- A_mat[j, i] <- 0.25 # half-sibs
  }
}
A_mat <- A_mat + diag(0.01, n_geno) # ensure PD
rownames(A_mat) <- colnames(A_mat) <- levels(dat$geno)

# Assemble
met_example <- list(
  dat = dat,
  G_mat = G_mat,
  A_mat = A_mat,
  M_geno = M_geno,
  n_geno = n_geno,
  n_env = n_env
)

# Save
usethis::use_data(met_example, overwrite = TRUE)
