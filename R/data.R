#' Example Multi-Environment Trial (MET) Dataset
#'
#' A synthetic multi-environment trial dataset suitable for demonstrating
#' mixed model formulations. Contains yield data for 10 genotypes evaluated
#' across 4 environments with 3 replicates per environment, plus spatial
#' layout, genomic markers, and multiple response types.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{`dat`}{A data frame with 120 rows and the following columns:
#'     \describe{
#'       \item{`geno`}{Factor: genotype identifier (G1--G10)}
#'       \item{`env`}{Factor: environment identifier (E1--E4)}
#'       \item{`rep`}{Factor: replicate (R1--R3)}
#'       \item{`block`}{Factor: incomplete block within rep (B1--B6)}
#'       \item{`row`}{Factor: row position in field layout}
#'       \item{`col`}{Factor: column position in field layout}
#'       \item{`yield`}{Numeric: grain yield (continuous, primary response)}
#'       \item{`x_cov`}{Numeric: a continuous covariate}
#'       \item{`bin_y`}{Integer: binary response (0/1)}
#'       \item{`count_y`}{Integer: count response (Poisson-like)}
#'     }
#'   }
#'   \item{`G_mat`}{10x10 genomic relationship matrix (G = MM'/2pq)}
#'   \item{`A_mat`}{10x10 pedigree-based relationship matrix}
#'   \item{`M_geno`}{10x50 marker matrix (50 SNP markers for 10 genotypes)}
#'   \item{`n_geno`}{Integer: number of genotypes (10)}
#'   \item{`n_env`}{Integer: number of environments (4)}
#' }
#'
#' @details
#' The dataset is generated with known variance components:
#' \itemize{
#'   \item Genotype variance: sigma_g^2 = 4
#'   \item Environment variance: sigma_e^2 = 9 (varies by environment)
#'   \item GxE variance: sigma_ge^2 = 2
#'   \item Residual variance: sigma^2 = 1
#' }
#'
#' The genomic relationship matrix is computed from 50 simulated SNP markers
#' using the VanRaden method 1. The pedigree matrix uses a simple half-sib
#' structure.
#'
#' @examples
#' data(met_example)
#' str(met_example$dat)
#' head(met_example$dat)
#'
#' # Simple model (small budget for example purposes)
#' \dontrun{
#' # live fit -- needs a backend (greta Python/TF, INLA, or brms/Stan)
#' fit <- flexybayes(
#'   fixed  = yield ~ env,
#'   random = ~ geno,
#'   data   = met_example$dat,
#'   n_samples = 100, warmup = 100, chains = 1, verbose = FALSE
#' )
#' summary(fit)
#' }
#'
#' @source Generated synthetically for package demonstration purposes.
"met_example"
