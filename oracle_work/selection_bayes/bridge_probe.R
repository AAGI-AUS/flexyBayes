suppressPackageStartupMessages(library(glmnet))
set.seed(3); n <- 50; p <- 6
A <- matrix(rnorm(n*p), n, p); beta <- c(2, -1, 0, 0, 0.5, 0); z <- A %*% beta + rnorm(n)
tau <- 0.7
# posterior mode of N(0,tau^2) prior with unit-variance whitened noise: (A'A + I/tau^2)^-1 A'z
closed <- solve(crossprod(A) + diag(p)/tau^2, crossprod(A, z))
lam_review <- 1/(n*tau^2); lam_v4 <- 1/tau^2
g1 <- as.numeric(coef(glmnet(A, z, alpha = 0, lambda = lam_review, standardize = FALSE, intercept = FALSE, thresh = 1e-14))[-1])
g2 <- as.numeric(coef(glmnet(A, z, alpha = 0, lambda = lam_v4, standardize = FALSE, intercept = FALSE, thresh = 1e-14))[-1])
cat("max|closed - glmnet(1/(n tau^2))| =", max(abs(closed - g1)), "\n")
cat("max|closed - glmnet(1/tau^2)|     =", max(abs(closed - g2)), "\n")
# lasso: Laplace scale b -> lambda = 1/(n b); check KKT of glmnet solution against objective 0.5||z-Ab||^2 + |b|_1 / b
b_scale <- 0.05; lam <- 1/(n*b_scale)
l1 <- as.numeric(coef(glmnet(A, z, alpha = 1, lambda = lam, standardize = FALSE, intercept = FALSE, thresh = 1e-14))[-1])
grad <- -crossprod(A, z - A %*% l1)  # gradient of 0.5||.||^2
kkt <- ifelse(l1 != 0, abs(grad + sign(l1)/b_scale), pmax(abs(grad) - 1/b_scale, 0))
cat("lasso KKT residual under 1/(n b) map, max =", max(kkt), "\n")
# glmnet default std/intercept ON silently changes the answer:
g3 <- as.numeric(coef(glmnet(A, z, alpha = 0, lambda = lam_review))[-1])
cat("max|closed - glmnet defaults| =", max(abs(closed - g3)), "\n")
cat("n =", n, " glmnet", as.character(packageVersion("glmnet")), "\n")
