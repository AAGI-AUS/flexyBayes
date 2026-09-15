suppressPackageStartupMessages({library(brms); library(projpred)})
options(brms.backend = "cmdstanr")
set.seed(2)
n <- 120; g <- factor(rep(1:12, each = 10))
X <- matrix(rnorm(n * 8), n, 8); colnames(X) <- paste0("x", 1:8)
u <- rnorm(12, 0, 0.7)
y <- 1 + 1.5 * X[,1] - 1 * X[,2] + u[g] + rnorm(n, 0, 0.5)
d <- data.frame(y, p1 = rnorm(n), X, g)
K <- diag(12); rownames(K) <- colnames(K) <- levels(g)
try_probe <- function(label, expr) {
  r <- tryCatch({v <- expr; list(ok = TRUE, v = v)}, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  cat(sprintf("\n=== %s: %s\n", label, if (r$ok) "OK" else paste("FAIL:", r$msg))); invisible(r)
}
t0 <- Sys.time()
m1 <- brm(y ~ p1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + (1|g), d,
          prior = prior(horseshoe(par_ratio = 0.3), class = "b"),
          chains = 2, iter = 1500, refresh = 0, seed = 1, silent = 2)
cat("fit m1:", format(difftime(Sys.time(), t0, units = "secs")), "\n")
r1 <- try_probe("P1 get_refmodel brmsfit (1|g) horseshoe", get_refmodel(m1))
if (r1$ok) {
  t0 <- Sys.time()
  vs <- try_probe("P2 varsel forward, nterms_max=6", varsel(r1$v, method = "forward", nterms_max = 6, verbose = FALSE, seed = 1))
  cat("varsel:", format(difftime(Sys.time(), t0, units = "secs")), "\n")
  if (vs$ok) { print(ranking(vs$v)); print(suggest_size(vs$v)); print(summary(vs$v, stats = "elpd")$perf_sub[, c("size","elpd","elpd.se")] ) }
  # grouped k-fold CV with custom folds (units = g), validated search
  folds <- as.integer(g) %% 4 + 1
  t0 <- Sys.time()
  cvs <- try_probe("P3 cv_varsel kfold, custom grouped folds, validate_search=TRUE",
    cv_varsel(r1$v, cv_method = "kfold", K = 4, folds = folds, method = "forward", nterms_max = 4,
              validate_search = TRUE, verbose = FALSE, seed = 1))
  cat("cv_varsel:", format(difftime(Sys.time(), t0, units = "secs")), "\n")
  if (cvs$ok) print(summary(cvs$v, stats = "elpd")$perf_sub[, c("size","elpd","elpd.se")])
  # search_terms to force-include p1 (protected) — every candidate model must contain p1
  st <- try_probe("P4 varsel with search_terms forcing p1", varsel(r1$v, method = "forward", nterms_max = 4, verbose = FALSE, seed = 1,
    search_terms = c("p1", paste0("p1 + ", c(paste0("x", 1:8), "(1 | g)")))))
  if (st$ok) print(ranking(st$v))
}
# gr(cov=) model
m2 <- brm(y ~ p1 + x1 + x2 + (1|gr(g, cov = K)), d, data2 = list(K = K), chains = 1, iter = 600, refresh = 0, seed = 1, silent = 2)
r2 <- try_probe("P5 get_refmodel with gr(g, cov = K)", get_refmodel(m2))
if (r2$ok) try_probe("P5b varsel with gr(cov)", varsel(r2$v, nterms_max = 2, verbose = FALSE, seed = 1))
# nlpar model
m3 <- brm(bf(y ~ prot + cand, prot ~ 0 + Intercept + p1 + (1|g), cand ~ 0 + x1 + x2 + x3, nl = TRUE), d,
          prior = c(prior(normal(0, 10), class = "b", nlpar = "prot"), prior(horseshoe(par_ratio = 0.3), class = "b", nlpar = "cand")),
          chains = 1, iter = 600, refresh = 0, seed = 1, silent = 2)
try_probe("P6 get_refmodel nlpar partition", get_refmodel(m3))
# log_lik / posterior_epred on new group
try_probe("P7 posterior_epred new group allow_new_levels", dim(posterior_epred(m1, newdata = transform(d[1:3,], g = factor("99")), allow_new_levels = TRUE, sample_new_levels = "gaussian")))
