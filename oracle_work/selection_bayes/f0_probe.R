suppressPackageStartupMessages(library(brms))
options(brms.backend = "cmdstanr")
set.seed(1)
n <- 60; d <- data.frame(y = rnorm(n), p1 = rnorm(n), p2 = rnorm(n),
  x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), g = factor(rep(1:6, 10)))
probe <- function(label, expr) {
  r <- tryCatch({ code <- expr; list(ok = TRUE, code = code) },
                error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  cat(sprintf("\n=== %s: %s\n", label, if (r$ok) "OK" else paste("FAIL:", r$msg)))
  if (r$ok) { l <- strsplit(r$code, "\n")[[1]]; cat(grep("hs|lambda|tau|b\\[|zb|scale_b|double_exp|R2D2|target \\+=|lprior", l, value = TRUE), sep = "\n") }
  invisible(r)
}
# 1 class-wide horseshoe with random intercept
probe("A class-wide horseshoe + (1|g)", make_stancode(y ~ p1 + x1 + x2 + x3 + (1|g), d,
  prior = prior(horseshoe(par_ratio = 0.3), class = "b")))
# 2 coefficient-specific laplace on candidates, normal on protected
probe("B coef-specific double_exponential + protected normal", make_stancode(y ~ p1 + x1 + x2 + x3 + (1|g), d,
  prior = c(prior(normal(0, 10), class = "b", coef = "p1"),
            prior(double_exponential(0, 0.5), class = "b", coef = "x1"),
            prior(double_exponential(0, 0.5), class = "b", coef = "x2"),
            prior(double_exponential(0, 0.5), class = "b", coef = "x3"))))
# 3 horseshoe class-wide + protected coef normal (expected FAIL)
probe("C horseshoe + coef override (expect FAIL)", make_stancode(y ~ p1 + x1 + x2 + x3, d,
  prior = c(prior(horseshoe(), class = "b"), prior(normal(0, 10), class = "b", coef = "p1"))))
# 4 nlpar partition
probe("D nlpar partition horseshoe/normal + (1|g)", make_stancode(
  bf(y ~ prot + cand, prot ~ 0 + Intercept + p1 + (1|g), cand ~ 0 + x1 + x2 + x3, nl = TRUE), d,
  prior = c(prior(normal(0, 10), class = "b", nlpar = "prot"),
            prior(horseshoe(par_ratio = 0.3), class = "b", nlpar = "cand"))))
# 5 R2D2 class-wide
probe("E class-wide R2D2", make_stancode(y ~ p1 + x1 + x2 + x3 + (1|g), d,
  prior = prior(R2D2(mean_R2 = 0.5, prec_R2 = 2), class = "b")))
# 6 elastic net via check=FALSE on indexed b, with stanvars index
sv <- stanvar(x = c(2L,3L,4L), name = "cand_idx", block = "data", scode = "array[3] int cand_idx;") +
      stanvar(scode = "real<lower=0> en_lambda = 0.5; real<lower=0,upper=1> en_alpha = 0.5;", block = "tdata")
probe("F elastic net target+= on b[cand_idx]", make_stancode(y ~ p1 + x1 + x2 + x3 + (1|g), d, stanvars = sv,
  prior = c(prior(normal(0, 10), class = "b", coef = "p1"),
            prior("target += -60 * en_lambda * (en_alpha * sum(abs(b[cand_idx])) + (1 - en_alpha) / 2 * dot_self(b[cand_idx]))", check = FALSE))))
# 7 does gr(cov=) coexist with horseshoe
K <- diag(6); rownames(K) <- colnames(K) <- levels(d$g)
probe("G horseshoe + gr(g, cov=K)", make_stancode(y ~ p1 + x1 + x2 + x3 + (1|gr(g, cov = K)), d, data2 = list(K = K),
  prior = prior(horseshoe(par_ratio = 0.3), class = "b")))
# 8 compile + sample smoke, timing
t0 <- Sys.time()
m <- brm(y ~ p1 + x1 + x2 + x3 + (1|g), d, prior = prior(horseshoe(par_ratio = 0.3), class = "b"),
         chains = 2, iter = 1000, refresh = 0, seed = 1, silent = 2)
cat("\n=== H smoke sample: ", format(difftime(Sys.time(), t0, units = "secs")), "\n")
print(max(rhat(m), na.rm = TRUE)); print(sum(nuts_params(m, pars = "divergent__")$Value))
print(class(m$fit)); print(cmdstanr::cmdstan_version())
