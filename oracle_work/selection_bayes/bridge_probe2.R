suppressPackageStartupMessages(library(glmnet))
set.seed(3); n <- 50; p <- 6
A <- matrix(rnorm(n*p), n, p); beta <- c(2, -1, 0, 0, 0.5, 0); z <- A %*% beta + rnorm(n)
tau <- 0.7; lam <- 1/(n*tau^2)
closed <- solve(crossprod(A) + diag(p)/tau^2, crossprod(A, z))
f_single <- glmnet(A, z, alpha = 0, lambda = lam, standardize = FALSE, intercept = FALSE, control = list(thresh = 1e-16, maxit = 1e6))
cat("single lambda, control thresh:", max(abs(closed - as.numeric(coef(f_single))[-1])), "\n")
path <- exp(seq(log(lam*100), log(lam/100), length.out = 100))
f_path <- glmnet(A, z, alpha = 0, lambda = path, standardize = FALSE, intercept = FALSE, control = list(thresh = 1e-16, maxit = 1e6))
cat("path + exact=TRUE:", max(abs(closed - as.numeric(coef(f_path, s = lam, exact = TRUE, x = A, y = z))[-1])), "\n")
# hypothesis: glmnet standardises y internally by sd(y) so lambda is on the y-sd scale for ridge? test lam * sd(z)?
sy <- sqrt(mean((z - mean(z))^2))
for (mult in c(1, sy, 1/sy)) {
  f <- glmnet(A, z, alpha = 0, lambda = lam * mult, standardize = FALSE, intercept = FALSE, control = list(thresh = 1e-16, maxit = 1e6))
  cat(sprintf("lambda x %.4f : %.3e\n", mult, max(abs(closed - as.numeric(coef(f))[-1]))))
}
