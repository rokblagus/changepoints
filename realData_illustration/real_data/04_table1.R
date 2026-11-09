##############################################################################
# Authors: Ana Radulovic and Rok Blagus
#
# 04_table1.R - Table 1 (main text) and the full table (Supplement S2)
#
# Combines results/sub/<method>/task*.csv (300 subsamples) and results/complete_*.csv.
# Output (folder results/):
#   table1_mean_sd.csv       mean (SD) of the selected K by method and subsample size,
#                            and the K selected on the complete subset (22 methods)
#   table1_distribution.csv  % of subsamples selecting K = 0, ..., 5
#   s_vs_s2_agreement.csv    % of subsamples where S and S2 select the same K
##############################################################################

rd <- function(what) {
  f <- list.files(file.path("results/sub", what), pattern = "^task[0-9]+\\.csv$", full.names = TRUE)
  if (length(f) == 0) return(NULL)
  do.call(rbind, lapply(f, read.csv))
}
long <- function(id, npy, K, m) data.frame(id = id, npy = npy, method = m, K = K)

L <- list()
nonp <- rd("nonperm")
if (!is.null(nonp)) {
  nm <- c(AIC = "AIC", BIC = "BIC", AIC_S = "AIC-S", BIC_S = "BIC-S",
          Davies_F = "Davies-F", Davies_R = "Davies-R")
  for (v in names(nm)) L[[length(L) + 1]] <- long(nonp$id, nonp$npy, nonp[[v]], nm[[v]])
}
label <- c(F_simulated = "F-S", F_simulated_type2 = "F-S2", F_permuted = "F-P",
           F_permuted_hetero = "F-P-het", gof_simulated = "GoF-S", gof_simulated_type2 = "GoF-S2",
           gof_permuted = "GoF-P", gof_permuted_hetero = "GoF-P-het")
for (t in names(label)) {
  r <- rd(t)
  if (is.null(r)) next
  L[[length(L) + 1]] <- long(r$id, r$npy, r$K,   label[[t]])
  L[[length(L) + 1]] <- long(r$id, r$npy, r$K_B, paste0(label[[t]], "-B"))
}
all <- do.call(rbind, L)
cat("Subsamples per method and size (should be 100 everywhere):\n")
print(table(all$method, all$npy))

order_m <- c("AIC", "BIC", "AIC-S", "BIC-S", "Davies-F", "Davies-R",
             "F-S", "F-S-B", "F-S2", "F-S2-B", "F-P", "F-P-B", "F-P-het", "F-P-het-B",
             "GoF-S", "GoF-S-B", "GoF-S2", "GoF-S2-B", "GoF-P", "GoF-P-B", "GoF-P-het", "GoF-P-het-B")
order_m <- order_m[order_m %in% all$method]

# selected K on the complete subset
complete <- setNames(rep(NA, length(order_m)), order_m)
if (file.exists("results/complete_nonperm.csv")) {
  cn <- read.csv("results/complete_nonperm.csv")
  complete[cn$method[cn$method %in% order_m]] <- cn$K[cn$method %in% order_m]
}
for (t in names(label)) {
  f <- sprintf("results/complete_%s_bw1.csv", t)
  if (!file.exists(f)) next
  x <- read.csv(f)
  if (label[[t]] %in% order_m) complete[label[[t]]] <- x$K
  if (paste0(label[[t]], "-B") %in% order_m) complete[paste0(label[[t]], "-B")] <- x$K_B
}

ms <- function(k) sprintf("%.1f (%.2f)", mean(k, na.rm = TRUE), sd(k, na.rm = TRUE))
t1 <- do.call(rbind, lapply(order_m, function(m) {
  d <- all[all$method == m, ]
  data.frame(method = m, n200 = ms(d$K[d$npy == 200]), n500 = ms(d$K[d$npy == 500]),
             n750 = ms(d$K[d$npy == 750]), complete = complete[[m]])
}))
write.csv(t1, "results/table1_mean_sd.csv", row.names = FALSE)

dist <- do.call(rbind, lapply(order_m, function(m) do.call(rbind, lapply(c(200, 500, 750), function(n) {
  k <- all$K[all$method == m & all$npy == n]
  data.frame(method = m, npy = n, t(sapply(0:5, function(j) round(100 * mean(k == j, na.rm = TRUE)))))
}))))
names(dist)[3:8] <- paste0("K", 0:5)
write.csv(dist, "results/table1_distribution.csv", row.names = FALSE)

# agreement between the two simulation-based variants
w <- reshape(all, idvar = c("id", "npy"), timevar = "method", direction = "wide")
agree <- function(a, b) if (all(paste0("K.", c(a, b)) %in% names(w)))
  round(100 * mean(w[[paste0("K.", a)]] == w[[paste0("K.", b)]], na.rm = TRUE), 1) else NA
ag <- data.frame(comparison = c("F-S vs F-S2", "F-S-B vs F-S2-B", "GoF-S vs GoF-S2", "GoF-S-B vs GoF-S2-B"),
                 percent_same_K = c(agree("F-S", "F-S2"), agree("F-S-B", "F-S2-B"),
                                    agree("GoF-S", "GoF-S2"), agree("GoF-S-B", "GoF-S2-B")))
write.csv(ag, "results/s_vs_s2_agreement.csv", row.names = FALSE)

cat("\nTable 1: mean (SD) of the selected number of changepoints\n"); print(t1, row.names = FALSE)
cat("\n"); print(ag, row.names = FALSE)
