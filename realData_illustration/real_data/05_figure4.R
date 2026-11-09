##############################################################################
# Authors: Ana Radulovic and Rok Blagus
# Part of the code for: "Estimating the number of changepoints in segmented
# regression: new sequential tests and comparative evaluation"
#
# 05_figure4.R - Figure 4: standardised annual mean BMI with the fitted models
#
# Grey points: annual means of BMI standardised for age and sex, using the
#   coefficients of a model with sex, age and calendar year as a categorical
#   variable (no assumption on the shape of the trend).
# Lines: predicted values for a child of mean age, weighted by the observed
#   sex proportions, from the segmented models with the numbers of changepoints
#   selected on the complete subset (locations from results/complete_models_by_K.csv,
#   i.e. as estimated by segmented in the selection procedure), and the model
#   selected by Davies-F (fixed locations from the sequential procedure).
# Requires: results of 02_complete_data.R. Output: results/figure4.pdf/.png
##############################################################################

library(ggplot2)
source("00_setup.R")
dat <- mydata$data
dat$girl <- as.numeric(dat$SEX == levels(dat$SEX)[2])

# methods shown in the main text
main <- c("AIC", "AIC-S", "BIC", "BIC-S", "Davies-R", "F-S", "F-P-het", "GoF-S", "GoF-P-het")
cn <- read.csv("results/complete_nonperm.csv")
selK <- setNames(cn$K, cn$method)
lab <- c(F_simulated = "F-S", F_permuted_hetero = "F-P-het",
         gof_simulated = "GoF-S", gof_permuted_hetero = "GoF-P-het")
for (t in names(lab)) {
  f <- sprintf("results/complete_%s_bw1.csv", t)
  if (file.exists(f)) selK[lab[[t]]] <- read.csv(f)$K
}
selK <- selK[names(selK) %in% main]
models <- read.csv("results/complete_models_by_K.csv")
parse_loc <- function(s) as.numeric(strsplit(as.character(s), ";")[[1]])

lines_def <- list()
for (k in sort(unique(selK[selK > 0]))) {
  lines_def[[sprintf("K = %d (%s)", k, paste(names(selK)[selK == k], collapse = ", "))]] <-
    parse_loc(models$locations[models$K == k])
}
dF <- cn[cn$method == "Davies-F", ]
if (dF$K > 0) lines_def[[sprintf("Davies-F (K = %d)", dF$K)]] <- sort(parse_loc(dF$locations))

# standardisation
m_cat  <- lm(Y ~ girl + AGE_YRS + factor(X), data = dat)
th_age <- coef(m_cat)[["AGE_YRS"]]; th_sex <- coef(m_cat)[["girl"]]
dat$Y_std <- dat$Y - th_age * (dat$AGE_YRS - mean(dat$AGE_YRS)) - th_sex * (dat$girl - mean(dat$girl))
means <- aggregate(Y_std ~ X, data = dat, FUN = mean)

# fitted lines (models with fixed locations; coefficients as in segmented at these locations)
xg <- seq(min(dat$X), max(dat$X), by = 0.05)
fl <- do.call(rbind, lapply(names(lines_def), function(nm) {
  tau <- lines_def[[nm]]
  U   <- matrix(sapply(tau, function(p) pmax(dat$X - p, 0)), ncol = length(tau))
  b   <- coef(lm(dat$Y ~ dat$girl + dat$AGE_YRS + dat$X + U))
  Ug  <- matrix(sapply(tau, function(p) pmax(xg - p, 0)), ncol = length(tau))
  data.frame(X = xg, model = nm,
             fitted = b[1] + b[2] * mean(dat$girl) + b[3] * mean(dat$AGE_YRS) + b[4] * xg +
                      as.vector(Ug %*% b[-(1:4)]))
}))
fl$model <- factor(fl$model, levels = names(lines_def))

p <- ggplot() +
  geom_point(data = means, aes(X, Y_std), colour = "grey55", size = 2.3) +
  geom_line(data = fl, aes(X, fitted, colour = model, linetype = model), linewidth = 1) +
  scale_linetype_manual(values = ifelse(grepl("^Davies-F", names(lines_def)), "dashed", "solid"),
                        name = NULL) +
  scale_colour_brewer(palette = "Set1", name = NULL) +
  scale_x_continuous(breaks = seq(1990, 2020, by = 5)) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", legend.direction = "vertical", panel.grid.minor = element_blank()) +
  labs(x = "Year", y = "BMI")
ggsave("results/figure4.pdf", p, width = 9, height = 6)
ggsave("results/figure4.png", p, width = 9, height = 6, dpi = 300)
write.csv(means, "results/figure4_annual_means.csv", row.names = FALSE)
print(lines_def)
cat("Saved: results/figure4.pdf, results/figure4.png\n")
