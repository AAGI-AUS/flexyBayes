# screen_repeatability.R -- can this structure return the same answer twice?
#
# Run before any further paired ASReml comparison. A single-run comparison is
# not evidence for a structure whose repeated fits disagree with each other:
# the number it produces is one draw from a distribution over modes, and
# comparing that draw against ASReml's stable answer measures the draw.
#
# Criterion, fixed here before the screen ran:
#
#   collapsed        a component's SD below 1% of the residual SD (the
#                    package's own .FB_COLLAPSE_FRACTION).
#   UNSTABLE         any component collapsed in at least one run and not in
#                    another, OR any non-collapsed component whose maximum
#                    across runs exceeds twice its minimum.
#   STABLE           neither condition fires.
#
# A factor of two is deliberately generous. INLA's integration is not
# deterministic and small movement is expected; a component that doubles
# between identical calls is not small movement.

.libPaths(path.expand(c(Sys.getenv("FB_LIB"),
                        "~/Library/R/arm64/4.5/library", .libPaths())))
suppressMessages(library(flexyBayes))
options(flexyBayes.silence_default_prior_note = TRUE)

N_RUNS <- as.integer(Sys.getenv("FB_N_RUNS", "3"))
COLLAPSE_FRACTION <- 0.01
INSTABILITY_RATIO <- 2

met <- function(n_env = NULL, n_gen = NULL) {
  e <- new.env()
  utils::data("besag.met", package = "agridat", envir = e)
  d <- e$besag.met
  d$env <- factor(d$county); d$gen <- factor(d$gen); d$rep <- factor(d$rep)
  d <- d[!is.na(d$yield), ]
  if (!is.null(n_env)) d <- d[d$env %in% levels(d$env)[seq_len(n_env)], ]
  if (!is.null(n_gen)) d <- d[d$gen %in% levels(d$gen)[seq_len(n_gen)], ]
  d$env <- droplevels(d$env); d$gen <- droplevels(d$gen)
  d <- d[order(d$env, d$gen), ]
  rownames(d) <- NULL
  d
}

sim_slope <- function() {
  set.seed(20260911L)
  s1 <- 240L
  g <- factor(rep(paste0("u", 1:24), each = 10L))
  x <- rep(seq(0, 9), times = 24L)
  b0 <- rnorm(24L, 0, 2)[as.integer(g)]
  b1 <- rnorm(24L, 0, 0.5)[as.integer(g)]
  data.frame(g = g, x = x, y = 10 + b0 + b1 * x + rnorm(s1, 0, 1))
}

structures <- list(
  "iid_intercept"   = list(d = function() met(), f = function(d)
    flexybayes(yield ~ env, random = ~ gen, data = d, backend = "inla", verbose = FALSE)),
  "nested"          = list(d = function() met(), f = function(d)
    flexybayes(yield ~ env, random = ~ env:rep, data = d, backend = "inla", verbose = FALSE)),
  "diag_gxe"        = list(d = function() met(), f = function(d)
    flexybayes(yield ~ env, random = ~ diag(env):gen, data = d, backend = "inla", verbose = FALSE)),
  "dsum_residual"   = list(d = function() met(), f = function(d)
    flexybayes(yield ~ env, random = ~ gen, residual = ~ dsum(~ units | env),
               data = d, backend = "inla", verbose = FALSE)),
  "combined_met"    = list(d = function() met(), f = function(d)
    flexybayes(yield ~ env, random = ~ diag(env):gen, residual = ~ dsum(~ units | env),
               data = d, backend = "inla", verbose = FALSE)),
  "us_dim2"         = list(d = function() met(n_env = 2L), f = function(d)
    flexybayes(yield ~ env, random = ~ us(env):gen, data = d, backend = "inla", verbose = FALSE)),
  "us_dim3"         = list(d = function() met(n_env = 3L), f = function(d)
    flexybayes(yield ~ env, random = ~ us(env):gen, data = d, backend = "inla", verbose = FALSE)),
  "fa1"             = list(d = function() met(n_gen = 30L), f = function(d)
    flexybayes(yield ~ env, random = ~ fa(env, 1):gen, data = d, backend = "inla", verbose = FALSE)),
  "cor_slope"       = list(d = function() sim_slope(), f = function(d)
    flexybayes(y ~ x + (x | g), data = d, backend = "inla", verbose = FALSE)),
  "uncor_slope"     = list(d = function() sim_slope(), f = function(d)
    flexybayes(y ~ x + (x || g), data = d, backend = "inla", verbose = FALSE))
)

id <- Sys.getenv("FB_STRUCT")
out <- Sys.getenv("FB_SCREEN_OUT")
st <- structures[[id]]
if (is.null(st)) stop("unknown structure: ", id, call. = FALSE)

d <- st$d()
runs <- list()
for (i in seq_len(N_RUNS)) {
  f <- tryCatch(suppressWarnings(st$f(d)), error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%s\tERROR\t%s\n", id, substr(conditionMessage(f), 1, 90)),
        file = out, append = TRUE)
    quit(save = "no", status = 0)
  }
  v <- summary(f)$varcomp
  runs[[i]] <- stats::setNames(v$estimate, v$component)
  last_fit <- f
}

comps <- Reduce(intersect, lapply(runs, names))
comps <- comps[!startsWith(comps, "cor_") & !startsWith(comps, "rho_")]
resid <- vapply(runs, function(r) {
  s <- r[names(r) == "sigma"]
  if (length(s)) unname(s[[1]]) else stats::median(r, na.rm = TRUE)
}, numeric(1))

flip <- FALSE; ratio_max <- 1
for (cm in setdiff(comps, "sigma")) {
  vals <- vapply(runs, function(r) unname(r[[cm]]), numeric(1))
  collapsed <- vals < COLLAPSE_FRACTION * resid
  if (any(collapsed) && !all(collapsed)) flip <- TRUE
  live <- vals[!collapsed]
  if (length(live) > 1 && min(live) > 0) {
    ratio_max <- max(ratio_max, max(live) / min(live))
  }
}
# Amendment, 2026-09-11: a third verdict, added after the first screen run
# returned STABLE for `diag_gxe` -- the structure already measured as
# bistable over 24 fits. All three screen runs had collapsed alike, so there
# was no flip and no live component to compare, and a consistently wrong
# answer scored as a repeatable one. Repeatability is not correctness, and a
# criterion that cannot tell them apart has to say so rather than pass.
#
# The restart probe settles what three repeats cannot: if a restart from a
# mode where the pinned components carry real variance finds a materially
# better marginal likelihood, the collapse is the optimiser's, not the data's.
all_collapsed <- FALSE
for (cm in setdiff(comps, "sigma")) {
  vals <- vapply(runs, function(r) unname(r[[cm]]), numeric(1))
  if (all(vals < COLLAPSE_FRACTION * resid)) all_collapsed <- TRUE
}
verdict <- if (flip || ratio_max > INSTABILITY_RATIO) {
  "UNSTABLE"
} else if (all_collapsed) {
  "ALL_COLLAPSED"
} else {
  "STABLE"
}
# The probe verdict from the last fit, which is the cheaper and sharper
# signal: one fit plus one restart, rather than three fits compared.
probe <- NULL
if (exists("last_fit", inherits = FALSE)) {
  probe <- last_fit$num_check$restart_probe
}
probe_txt <- if (is.null(probe)) {
  "probe=not_run"
} else {
  sprintf("probe=%s(delta=%+.1f)", if (isTRUE(probe$better)) "FLAGGED" else "clear", probe$delta)
}
cat(sprintf("%s\t%s\tn_comp=%d\tmax_ratio=%.2f\tcollapse_flip=%s\t%s\truns=%d\n",
            id, verdict, length(comps), ratio_max, flip, probe_txt, N_RUNS),
    file = out, append = TRUE)
cat(sprintf("%-16s %-14s max_ratio=%6.2f flip=%-5s %s\n",
            id, verdict, ratio_max, flip, probe_txt))
