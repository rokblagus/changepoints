##############################################################################
# Authors: Ana Radulovic and Rok Blagus
# Part of the code for: "Estimating the number of changepoints in segmented
# regression: new sequential tests and comparative evaluation"
#
# 06_heteroscedasticity.R - does the residual variability depend on the year?
#                           (justification of the heteroscedasticity-adjusted tests)
#
# Residuals from a model with sex, age and calendar year as a categorical
# variable (no assumption on the shape of the trend).
#  1. SD of the residuals by year with 95% CIs, and the Nadaraya-Watson
#     estimate of the SD used by the hetero-adjusted tests (bandwidth c/sqrt(n),
#     c = (range of year)^2 / 2, as in the analysis);
#  2. SD of the residuals by deciles of the fitted values;
#  3. Breusch-Pagan (Koenker) test and Brown-Forsythe test across years.
# Output (folder results/): hetero_by_year.csv, hetero_tests.csv,
#   hetero_by_year.pdf/.png, hetero_by_fitted.pdf/.png
##############################################################################

library(ggplot2)
source("00_setup.R")
d <- mydata$data
n <- nrow(d)

m <- lm(Y ~ SEX + AGE_YRS + factor(X), data = d)
d$r <- resid(m); d$fit <- fitted(m)

sd_ci <- function(r) {
  k <- length(r); s <- sd(r)
  c(n = k, sd = s, lower = s * sqrt((k - 1) / qchisq(0.975, k - 1)),
    upper = s * sqrt((k - 1) / qchisq(0.025, k - 1)))
}
by_year <- do.call(rbind, lapply(split(d$r, d$X), sd_ci))
by_year <- data.frame(X = as.numeric(rownames(by_year)), by_year)
write.csv(by_year, "results/hetero_by_year.csv", row.names = FALSE)

h <- (diff(range(d$X))^2 / 2) / sqrt(n)
d$sd_kernel <- sqrt(sigma_u_by_x(d$X, d$r, h))
kern <- unique(d[, c("X", "sd_kernel")]); kern <- kern[order(kern$X), ]

p1 <- ggplot(by_year, aes(X, sd)) +
  geom_hline(yintercept = sd(d$r), linetype = "dotted", colour = "grey40") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, colour = "grey55") +
  geom_point(size = 2) +
  geom_line(data = kern, aes(X, sd_kernel), colour = "#E41A1C", linewidth = 1) +
  scale_x_continuous(breaks = seq(1990, 2020, by = 5)) +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank()) +
  labs(x = "Year", y = "SD of residuals")
ggsave("results/hetero_by_year.pdf", p1, width = 8, height = 5)
ggsave("results/hetero_by_year.png", p1, width = 8, height = 5, dpi = 300)

d$dec <- cut(d$fit, quantile(d$fit, 0:10 / 10), include.lowest = TRUE)
by_fit <- do.call(rbind, lapply(split(d, d$dec), function(z) c(fit = mean(z$fit), sd_ci(z$r))))
by_fit <- as.data.frame(by_fit)
p2 <- ggplot(by_fit, aes(fit, sd)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.05, colour = "grey55") +
  geom_point(size = 2) + geom_line() +
  theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank()) +
  labs(x = "Fitted BMI (decile means)", y = "SD of residuals")
ggsave("results/hetero_by_fitted.pdf", p2, width = 8, height = 5)
ggsave("results/hetero_by_fitted.png", p2, width = 8, height = 5, dpi = 300)

bp <- n * summary(lm(I(r^2) ~ factor(X), data = d))$r.squared
df_bp <- length(unique(d$X)) - 1
d$absdev <- abs(d$r - ave(d$r, d$X, FUN = median))
bf <- anova(lm(absdev ~ factor(X), data = d))
tests <- data.frame(
  test      = c("Breusch-Pagan (Koenker), year as factor", "Brown-Forsythe, by year"),
  statistic = c(round(bp, 1), round(bf$`F value`[1], 2)),
  df        = c(df_bp, paste(bf$Df[1], bf$Df[2], sep = ", ")),
  p_value   = c(signif(pchisq(bp, df_bp, lower.tail = FALSE), 3), signif(bf$`Pr(>F)`[1], 3)))
write.csv(tests, "results/hetero_tests.csv", row.names = FALSE)

print(round(by_year, 3), row.names = FALSE)
cat(sprintf("\nSmallest SD: %.2f (%d) | largest SD: %.2f (%d)\n",
            min(by_year$sd), by_year$X[which.min(by_year$sd)],
            max(by_year$sd), by_year$X[which.max(by_year$sd)]))
print(tests, row.names = FALSE)
