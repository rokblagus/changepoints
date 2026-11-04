
##this is for explaining the conservative nature of P approaches.


rm(list = ls())

# detach("package:segmented", unload = TRUE)
# detach("package:MASS", unload = TRUE)

 
delta0 <- c(0.2)
# delta1 <- seq(from = 0.2, to = 1, by = 0.4) # issues with floats!
delta1 <- c(0.6)
task <- c("medium")

# var.sigma <- c(0.5, 1, 1.5)
var.sigma <- c(3)

model <- "linear"
phi.beta <- NULL
intercept <- 1
mean_y <- NULL

start <- 0
no <- 2:3
n <- c(250, 500)

pdf("figures/plot_conser.pdf",height=9,width = 9)
par(mfrow = c(2, 2))

for (i in delta0) {
  for (ii in delta1) {
    for (j in task) {
      for (jj in no) {

        if (jj == 1) alteri <- FALSE else alteri <- TRUE

        for (z in n) {
          for (zz in alteri) {
            for (mm in var.sigma) {

              nm <- paste(
                "d0", i, "d1", ii, j,
                "a", ifelse(zz == TRUE, 1, 0),
                "no", jj, "n", z, model,
                sep = "."
              )

              if (model == "linear" | model == "quantile") {
                nm <- paste(nm, "vs", mm, "in", intercept, sep = ".")
              } else {
                nm <- paste(nm, "mY", mean_y, sep = ".")
              }

              if (model == "beta") {
                nm <- paste(nm, "ph", phi.beta, sep = ".")
              }

              nm <- paste(nm, ".txt", sep = "")
              nmr1 <- paste("results/r1", nm, sep = ".")

              start <- start + 1

              dd1 <- read.table(
                nmr1,
                sep = "\t",
                fill = TRUE,
                header = FALSE
              )

              names(dd1) <- c(
                "F-P-set-p", "F-P-set-ptilda", "F-P-set-ps", "F-P-set-F",
                "F-P-p", "F-P-ptilda", "F-P-ps", "F-P-F",
                "GoF-P-p", "GoF-P-ptilda", "GoF-P-ps", "GoF-P-F",
                "GoF-P-set-p", "GoF-P-set-ptilda", "GoF-P-set-ps", "GoF-P-set-F"
              )

              # --------------------------------------------------
              # Boxplots with all observations shown as jittered
              # light-gray points
              # --------------------------------------------------

              dd_min <- (dd1[, c(4, 8, 12, 16)] + 1) / 1001

              names(dd_min) <- c(
                "F-P-set", "F-P", "GoF-P", "GoF-P-set"
              )

              # Preserve the original plotting order
              dd_min <- dd_min[, c(1, 2, 4, 3)]

              # Calculate a valid numeric y-axis maximum
              ymax <- max(0.06, unlist(dd_min), na.rm = TRUE)

              # Draw boxplots first to establish axes
              boxplot(
                dd_min,
                main = paste0("n = ", z, ", K0 = ", jj),
                ylab = "(F + 1)/(B + 1)",
                outline = FALSE,
                ylim = c(0, ymax),
                col = "white"
              )

              # Add individual observations with horizontal jitter
              stripchart(
                dd_min,
                vertical = TRUE,
                method = "jitter",
                jitter = 0.18,
                pch = 16,
                cex = 0.65,
                col = grDevices::adjustcolor(
                  "gray50",
                  alpha.f = 0.4
                ),
                add = TRUE
              )

              # Overlay boxplots on the points
              boxplot(
                dd_min,
                add = TRUE,
                outline = FALSE,
                axes = FALSE,
                col = grDevices::adjustcolor(
                  "white",
                  alpha.f = 0.75
                ),
                border = "black"
              )

              # Nominal significance level
              abline(h = 0.05, lty = "dashed")

              # --------------------------------------------------
              # Summary of empirical rejection rates
              # --------------------------------------------------

              my.sum <- function(x) {
                res <- c(
                  length(x),
                  mean(x <= 0.05, na.rm = TRUE),
                  sum(is.na(x))
                )

                names(res) <- c("Nsim", "rej_05", "nNA")
                res
              }

              res <- t(
                apply(
                  dd1[, -c(4, 8, 12, 16)],
                  2,
                  my.sum
                )
              )

              # --------------------------------------------------
              # Summary of the fourth statistic for each method
              # --------------------------------------------------

              my.sum2 <- function(x) {
                res <- c(
                  mean(x, na.rm = TRUE),
                  min(x, na.rm = TRUE),
                  max(x, na.rm = TRUE),
                  sum(is.na(x)),
                  length(x)
                )

                names(res) <- c(
                  "meanF", "minF", "maxF", "nNA", "Nsim"
                )

                res
              }

              res2 <- t(
                apply(
                  dd1[, c(4, 8, 12, 16)],
                  2,
                  my.sum2
                )
              )

              # --------------------------------------------------
              # Store summaries and simulation settings
              # --------------------------------------------------

              dfi1 <- as.data.frame(res)
              dfi1$method <- rownames(dfi1)

              dfi1$delta0 <- i
              dfi1$delta1 <- ii
              dfi1$task <- j
              dfi1$no <- jj
              dfi1$n <- z
              dfi1$alter <- zz
              dfi1$var.sigma <- mm

              dfi12 <- as.data.frame(res2)
              dfi12$method <- rownames(dfi12)

              dfi12$delta0 <- i
              dfi12$delta1 <- ii
              dfi12$task <- j
              dfi12$no <- jj
              dfi12$n <- z
              dfi12$alter <- zz
              dfi12$var.sigma <- mm

              if (start == 1) {
                df <- dfi1
                df2 <- dfi12
              } else {
                df <- rbind(df, dfi1)
                df2 <- rbind(df2, dfi12)
              }

            }
          }
        }
      }
    }
  }
}

dev.off()



library(xtable)

# --------------------------------------------------
# Prepare empirical rejection rates, MCSE, and NA percentages
# --------------------------------------------------

df_table <- df

# Number of valid simulation replicates
df_table$Nvalid <- df_table$Nsim - df_table$nNA

# MCSE on the probability scale
df_table$MCSE <- with(
  df_table,
  sqrt(rej_05 * (1 - rej_05) / Nvalid)
)

# Express rejection rates and MCSEs as percentages
df_table$rate_pct <- 100 * df_table$rej_05
df_table$MCSE_pct <- 100 * df_table$MCSE

# Percentage of simulation runs with missing p-values
df_table$na_pct <- 100 * df_table$nNA / df_table$Nsim

# Format rejection rate and MCSE on the first line,
# and NA percentage on the second line
df_table$rate_mcse <- ifelse(
  is.na(df_table$rate_pct) | is.na(df_table$MCSE_pct),
  "NA",
  paste0(
    "\\parbox[t]{1.8cm}{\\centering ",
    sprintf("%.2f (%.2f)", df_table$rate_pct, df_table$MCSE_pct),
    "\\\\[-3pt]{\\tiny NA: ",
    sprintf("%.1f", df_table$na_pct),
    "}\\\\[3pt]}"
  )
)

# Identify each simulation setting
df_table$setting <- paste0(
  "n=", df_table$n,
  ", K0=", df_table$no
)

# Keep the columns needed for the table
table_wide <- reshape(
  df_table[, c("method", "setting", "rate_mcse")],
  idvar = "method",
  timevar = "setting",
  direction = "wide"
)

# Clean column names
names(table_wide) <- sub(
  "^rate_mcse\\.",
  "",
  names(table_wide)
)

names(table_wide)[1] <- "Method"

# --------------------------------------------------
# Create LaTeX table
# --------------------------------------------------

align_vec <- c(
  "l",  # row-name alignment required by xtable
  "l",  # Method
  rep("p{1.8cm}", ncol(table_wide) - 1)
)

tab <- xtable(
  table_wide,
  caption = paste(
    "Empirical rejection rates (\\%) at the 5\\% significance level.",
    "Monte Carlo standard errors, expressed in percentage points,",
    "are given in parentheses. Percentages of missing p-values are",
    "shown on the second line."
  ),
  label = "tab:empirical_rejection_rates",
  align = align_vec
)

latex_tab <- print(
  tab,
  include.rownames = FALSE,
  sanitize.text.function = identity,
  booktabs = TRUE,
  comment = FALSE,
  floating = TRUE,
  print.results = FALSE
)

# Reduce horizontal spacing between columns for this table only
latex_tab <- sub(
  "\\\\begin\\{tabular\\}",
  "\\\\setlength{\\\\tabcolsep}{3pt}\n\\\\begin{tabular}",
  latex_tab
)

cat(latex_tab)


#see where we are

plot(df$Nsim)



