##############################################################################
# Authors: Ana Radulovic and Rok Blagus
# Part of the code for: "Estimating the number of changepoints in segmented
# regression: new sequential tests and comparative evaluation"
#
# 00_setup.R - common setup for all scripts of the real-data illustration
#
#  - loads the anonymised subset (created by 01_make_subset.R) in the format
#    used by the functions: mydata$data with columns Y (BMI), SEX, AGE_YRS
#    (covariates) and X (calendar year, the segmented predictor);
#  - defines the sequential testing procedure, the subsampling design for
#    Table 1, and the list of all jobs run by the workers (02_worker.R).
##############################################################################

for (pkg in c("segmented", "data.table"))
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
setDTthreads(1)                 # one core per R process

MAX_SEG   <- 5                  # K_max, as in the simulation study
ALPHA     <- 0.05               # significance level (ALPHA / MAX_SEG with Bonferroni)
B         <- 1000               # number of permutations / simulations per p-value
DATA_FILE <- "slofit_subset.Rdata"
for (d in c("results", "results/sub", "logs", "locks")) dir.create(d, showWarnings = FALSE)

cat("R:", R.version.string, "| segmented:", as.character(packageVersion("segmented")),
    "| RNG:", paste(RNGkind(), collapse = "/"), "\n")
if (packageVersion("segmented") != "2.2.2")
  warning("The results in the paper were obtained with segmented 2.2-2")

if (!file.exists(DATA_FILE)) stop("Data file not found: run 01_make_subset.R first")
load(DATA_FILE)                 # object 'slofit': BMI (with added noise), SEX, AGE_YRS, YEAR
mydata <- list(data = data.frame(Y = slofit$BMI, SEX = slofit$SEX,
                                 AGE_YRS = slofit$AGE_YRS, X = slofit$YEAR))

source("functions/functions_ic_davies.R")   # AIC, BIC, AIC-S, BIC-S, Davies-F, Davies-R
source("functions/functions_tests.R")       # F-type and GoF tests

# The eight tests (each is also evaluated with the Bonferroni-adjusted level)
TESTS <- c("F_simulated", "F_simulated_type2", "F_permuted", "F_permuted_hetero",
           "gof_simulated", "gof_simulated_type2", "gof_permuted", "gof_permuted_hetero")

# ---------------------------------------------------------------------------
# Sequential procedure: p-values for H0: K = 0, 1, ..., K_max - 1.
# A seed is set before every step. stop_early = TRUE stops at the first
# p > ALPHA (used for the subsamples: this gives the same selected K, also with
# the Bonferroni adjustment, and saves time); on the complete subset all five
# p-values are computed.
# ---------------------------------------------------------------------------
seq_test <- function(test, dat, stop_early, verbose = FALSE) {
  fun <- get(test)
  p <- rep(NA_real_, MAX_SEG)
  for (k in 0:(MAX_SEG - 1)) {
    set.seed(1000 + k)
    p[k + 1] <- tryCatch(fun(B, dat, k), error = function(e) NA)
    if (verbose) cat(sprintf("   H0: K = %d: p = %s\n", k, signif(p[k + 1], 4)))
    if (stop_early && (is.na(p[k + 1]) || p[k + 1] > ALPHA)) break
  }
  p
}

# Selected number of changepoints: the number of rejections before the first
# non-rejection; a missing p-value stops the procedure (as in the simulation study)
k_hat <- function(pv, alpha) {
  if (is.na(pv[1])) return(NA)
  k <- 0
  for (j in seq_len(MAX_SEG)) {
    if (is.na(pv[j]) || pv[j] > alpha) break
    k <- j
  }
  k
}
pstr <- function(p) paste(signif(p, 4), collapse = "; ")
loc  <- function(x) if (length(x) == 0) "" else paste(round(sort(x), 2), collapse = "; ")

# ---------------------------------------------------------------------------
# Subsampling design for Table 1: 100 replicates for each of 200, 500 and 750
# records per calendar year (300 tasks). Every task has its own seed, so the
# subsample does not depend on which process or computer draws it.
# ---------------------------------------------------------------------------
make_tasks <- function() {
  tasks <- expand.grid(rep = 1:100, npy = c(200, 500, 750))
  tasks$id <- seq_len(nrow(tasks))
  tasks
}
draw_subsample <- function(full_data, npy, rep) {
  idx_by_year <- split(seq_len(nrow(full_data)), full_data$X)
  set.seed(100000 * npy + rep)
  idx <- unlist(lapply(idx_by_year, function(ii) if (length(ii) <= npy) ii else sample(ii, npy)))
  full_data[sort(idx), ]
}

# ---------------------------------------------------------------------------
# All jobs of the analysis. Each job is one method (or the group of the six
# non-test methods, "nonperm") on the complete subset or on one subsample.
# The approximate cost (hours, from earlier runs) is used only to start the
# longest jobs first, which balances the load over the workers.
# ---------------------------------------------------------------------------
job_list <- function() {
  perm_cost <- rbind(F_permuted          = c(0.6, 3.8, 5.1),
                     F_permuted_hetero   = c(1.3, 3.9, 5.6),
                     gof_permuted        = c(0.5, 1.0, 1.8),
                     gof_permuted_hetero = c(0.8, 1.8, 2.6))
  colnames(perm_cost) <- c("200", "500", "750")
  J <- data.frame(type = "complete", what = c("nonperm", TESTS), bwf = 1,
                  id = NA, npy = NA, rep = NA,
                  cost = c(1, 0.2, 0.2, 12, 12, 0.5, 0.5, 8, 8))
  J <- rbind(J, data.frame(type = "complete",
                           what = rep(c("F_permuted_hetero", "gof_permuted_hetero"), 2),
                           bwf = rep(c(0.5, 2), each = 2), id = NA, npy = NA, rep = NA,
                           cost = c(12, 8, 12, 8)))
  tasks <- make_tasks()
  for (w in c("nonperm", TESTS)) {
    cost <- if (w %in% rownames(perm_cost)) perm_cost[w, as.character(tasks$npy)]
            else if (w == "nonperm") 0.3 * tasks$npy / 750 else 0.1 * tasks$npy / 750
    J <- rbind(J, data.frame(type = "sub", what = w, bwf = 1, id = tasks$id,
                             npy = tasks$npy, rep = tasks$rep, cost = as.numeric(cost)))
  }
  J$name <- ifelse(J$type == "complete", sprintf("complete_%s_bw%s", J$what, J$bwf),
                   sprintf("sub_%s_task%03d", J$what, J$id))
  J$name[J$type == "complete" & J$what == "nonperm"] <- "complete_nonperm"
  J[order(-J$cost, J$type, J$what, J$id), ]
}
result_file <- function(job) {
  if (job$type == "complete") file.path("results", paste0(job$name, ".csv"))
  else file.path("results/sub", job$what, sprintf("task%03d.csv", job$id))
}

# ---------------------------------------------------------------------------
# Runs one job and saves its result (one small file per job, so that many
# processes can run in parallel without writing to the same file)
# ---------------------------------------------------------------------------
run_job <- function(job) {
  t0   <- Sys.time()
  secs <- function() round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  out  <- result_file(job)
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  dat  <- if (job$type == "complete") mydata$data
          else draw_subsample(mydata$data, job$npy, job$rep)

  if (job$what == "nonperm") {
    sub <- list(data = dat)
    set.seed(1); ic <- tryCatch(ic_rok(sub, MAX_SEG), error = function(e) NULL)
    set.seed(1); dF <- tryCatch(davies_fixed(sub, MAX_SEG), error = function(e) NULL)
    set.seed(1); dR <- tryCatch(davies_reestimated(sub, MAX_SEG), error = function(e) NULL)
    g <- function(x, el) if (is.null(x)) NA else x[[el]]
    if (job$type == "complete") {
      # selected K with the changepoint locations, and all models with K = 0..5
      res <- data.frame(
        method    = c("AIC", "BIC", "AIC-S", "BIC-S", "Davies-F", "Davies-R"),
        K         = c(g(ic, "n.aic"), g(ic, "n.bic"), g(ic, "sel.aic"), g(ic, "sel.bic"),
                      g(dF, "nm.davies"), g(dR, "nm.davies")),
        locations = c(loc(ic$psi_list[[ic$n.aic + 1]]), loc(ic$psi_list[[ic$n.bic + 1]]),
                      loc(ic$psi_selaic), loc(ic$psi_selbic), loc(dF$chpnts), loc(dR$chpnts)))
      models <- data.frame(K = 0:MAX_SEG, locations = sapply(ic$psi_list, loc),
                           AIC = round(ic$aics, 2), BIC = round(ic$bics, 2))
      write.csv(models, "results/complete_models_by_K.csv", row.names = FALSE)
    } else {
      res <- data.frame(id = job$id, npy = job$npy, rep = job$rep, n = nrow(dat),
                        AIC = g(ic, "n.aic"), BIC = g(ic, "n.bic"),
                        AIC_S = g(ic, "sel.aic"), BIC_S = g(ic, "sel.bic"),
                        Davies_F = g(dF, "nm.davies"), Davies_R = g(dR, "nm.davies"))
    }
  } else {
    BW_FACTOR <<- job$bwf          # bandwidth factor for the hetero-adjusted tests
    p <- seq_test(job$what, dat, stop_early = (job$type == "sub"))
    BW_FACTOR <<- 1
    res <- data.frame(id = job$id, npy = job$npy, rep = job$rep, n = nrow(dat),
                      method = job$what, bw_factor = job$bwf,
                      K = k_hat(p, ALPHA), K_B = k_hat(p, ALPHA / MAX_SEG), p_values = pstr(p))
  }
  res$seconds <- secs()
  write.csv(res, out, row.names = FALSE)
  res
}
