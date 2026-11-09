##############################################################################
# Authors: Ana Radulovic and Rok Blagus
# Part of the code for: "Estimating the number of changepoints in segmented
# regression: new sequential tests and comparative evaluation"
#
# 07_coefficients.R - coefficients with 95% CIs for the selected models, and a
#                     check of the stability of the changepoint locations
#
# 1. For each number of changepoints selected on the complete subset (K > 0),
#    the segmented model is refitted from the locations estimated in the
#    selection procedure (results/complete_models_by_K.csv), so that the
#    coefficients correspond to the reported locations. Wald 95% CIs for the
#    regression coefficients and slope changes, confint() for the locations.
#    For Davies-F, the model with the fixed locations is fitted by lm(); its
#    CIs are conditional on these locations.
# 2. Local optima check: the segmented model with K changepoints is fitted
#    with 20 different seeds (bootstrap restarting, n.boot = 50), and the
#    distinct solutions are reported with their residual sum of squares.
# Output (folder results/): coefficients.csv, local_optima_check.csv
##############################################################################

source("00_setup.R")
dat  <- mydata$data
fit0 <- lm(Y ~ SEX + AGE_YRS + X, data = dat)

cn <- read.csv("results/complete_nonperm.csv")
selK <- cn$K
for (f in list.files("results", pattern = "^complete_.*_bw1\\.csv$", full.names = TRUE))
  selK <- c(selK, read.csv(f)$K)
Ks <- sort(unique(selK[!is.na(selK) & selK > 0]))
models <- read.csv("results/complete_models_by_K.csv")
parse_loc <- function(s) as.numeric(strsplit(as.character(s), ";")[[1]])

coef_table <- function(est, se, label, K) {
  data.frame(model = label, K = K, term = names(est), est = round(est, 4),
             lower = round(est - 1.96 * se, 4), upper = round(est + 1.96 * se, 4))
}

out <- NULL
for (k in Ks) {
  start <- parse_loc(models$locations[models$K == k])
  f <- segmented(fit0, seg.Z = ~X, psi = list(X = start),
                 control = seg.control(n.boot = 0, alpha = 0.01))
  est <- coef(f); se <- sqrt(diag(vcov(f)))[names(est)]
  ct  <- coef_table(est, se, sprintf("segmented, K = %d", k), k)
  ct  <- ct[!grepl("^psi", ct$term), ]
  ci  <- confint(f)
  ps  <- data.frame(model = sprintf("segmented, K = %d", k), K = k,
                    term = paste0("tau", seq_len(nrow(ci))), est = round(ci[, 1], 2),
                    lower = round(ci[, 2], 2), upper = round(ci[, 3], 2))
  out <- rbind(out, ct, ps)
}

dF <- cn[cn$method == "Davies-F", ]
if (!is.na(dF$K) && dF$K > 0) {
  tau <- sort(parse_loc(dF$locations))
  U <- matrix(sapply(tau, function(p) pmax(dat$X - p, 0)), ncol = length(tau)); colnames(U) <- paste0("U", seq_along(tau), ".X")
  fl <- lm(dat$Y ~ dat$SEX + dat$AGE_YRS + dat$X + U)
  est <- coef(fl); se <- sqrt(diag(vcov(fl)))
  out <- rbind(out, coef_table(est, se, sprintf("Davies-F, K = %d (fixed locations)", dF$K), dF$K),
               data.frame(model = sprintf("Davies-F, K = %d (fixed locations)", dF$K), K = dF$K,
                          term = paste0("tau", seq_along(tau)), est = round(tau, 2), lower = NA, upper = NA))
}
rownames(out) <- NULL
write.csv(out, "results/coefficients.csv", row.names = FALSE)
print(out, row.names = FALSE)

# local optima check
chk <- NULL
for (k in Ks) for (s in 1:20) {
  set.seed(s)
  f <- try(segmented(fit0, seg.Z = ~X, npsi = k, control = seg.control(n.boot = 50, alpha = 0.01)),
           silent = TRUE)
  if (inherits(f, "try-error") || is.null(f$psi) || nrow(f$psi) != k) next
  chk <- rbind(chk, data.frame(K = k, seed = s,
                               locations = paste(round(sort(f$psi[, "Est."]), 1), collapse = "; "),
                               RSS = round(sum(resid(f)^2), 2)))
}
sol <- aggregate(seed ~ K + locations + RSS, data = chk, FUN = length)
names(sol)[names(sol) == "seed"] <- "n_seeds"
sol <- sol[order(sol$K, sol$RSS), ]
write.csv(sol, "results/local_optima_check.csv", row.names = FALSE)
cat("\nDistinct solutions over 20 seeds:\n"); print(sol, row.names = FALSE)
