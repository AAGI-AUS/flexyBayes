suppressPackageStartupMessages({library(brms); library(glmnet); library(cmdstanr)})
set.seed(3); n <- 50; p <- 6
A <- matrix(rnorm(n*p), n, p); colnames(A) <- paste0("x", 1:p); beta <- c(2, -1, 0, 0, 0.5, 0); z <- as.numeric(A %*% beta + rnorm(n))
d <- data.frame(y = z, A)
b_scale <- 0.05
pri <- c(prior(constant(1), class = "sigma"), set_prior(sprintf("double_exponential(0, %g)", b_scale), class = "b"))
f <- bf(y ~ 0 + x1 + x2 + x3 + x4 + x5 + x6)
code <- make_stancode(f, d, prior = pri); sd <- make_standata(f, d, prior = pri)
mod <- cmdstan_model(write_stan_file(code), quiet = TRUE)
opt <- mod$optimize(data = sd, algorithm = "lbfgs", tol_rel_grad = 1e-12, tol_grad = 1e-12, tol_param = 1e-12, iter = 20000, seed = 1, refresh = 0)
b_map <- as.numeric(opt$mle(variables = "b"))
lam <- 1/(n*b_scale)
b_gl <- as.numeric(coef(glmnet(A, z, alpha = 1, lambda = lam, standardize = FALSE, intercept = FALSE, control = list(thresh = 1e-16, maxit = 1e6)))[-1])
cat("lasso MAP vs glmnet 1/(n b): max abs diff", max(abs(b_map - b_gl)), "\n"); print(rbind(map = b_map, glmnet = b_gl))
# ridge
tau <- 0.7
pri2 <- c(prior(constant(1), class = "sigma"), set_prior(sprintf("normal(0, %g)", tau), class = "b"))
code2 <- make_stancode(f, d, prior = pri2); sd2 <- make_standata(f, d, prior = pri2)
opt2 <- cmdstan_model(write_stan_file(code2), quiet = TRUE)$optimize(data = sd2, algorithm = "lbfgs", tol_rel_grad = 1e-12, tol_grad = 1e-12, tol_param = 1e-12, iter = 20000, seed = 1, refresh = 0)
closed <- as.numeric(solve(crossprod(A) + diag(p)/tau^2, crossprod(A, z)))
cat("ridge MAP vs closed form: max abs diff", max(abs(as.numeric(opt2$mle(variables = "b")) - closed)), "\n")
# elastic net via check=FALSE, sampled + MAP
alpha <- 0.5; lam_en <- 0.3
pri3 <- c(prior(constant(1), class = "sigma"),
          set_prior(sprintf("target += -%d * %g * (%g * sum(abs(b)) + %g / 2 * dot_self(b))", n, lam_en, alpha, 1 - alpha), check = FALSE))
code3 <- make_stancode(f, d, prior = pri3); sd3 <- make_standata(f, d, prior = pri3)
m3 <- cmdstan_model(write_stan_file(code3), quiet = TRUE)
opt3 <- m3$optimize(data = sd3, algorithm = "lbfgs", tol_rel_grad = 1e-12, tol_grad = 1e-12, tol_param = 1e-12, iter = 20000, seed = 1, refresh = 0)
b_gl3 <- as.numeric(coef(glmnet(A, z, alpha = alpha, lambda = lam_en, standardize = FALSE, intercept = FALSE, control = list(thresh = 1e-16, maxit = 1e6)))[-1])
cat("elastic net MAP vs glmnet(alpha, lambda) [direct objective match]: max abs diff", max(abs(as.numeric(opt3$mle(variables = "b")) - b_gl3)), "\n")
s3 <- m3$sample(data = sd3, chains = 2, iter_sampling = 1000, iter_warmup = 1000, seed = 1, refresh = 0, show_messages = FALSE)
print(s3$diagnostic_summary()); print(s3$summary("b")[, c("variable","mean","rhat","ess_bulk")])
