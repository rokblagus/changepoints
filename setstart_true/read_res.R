

delta0=c(0.2)
#delta1=seq(from=0.2,to=1,by=0.4) #issues with floats!!
delta1<-c(0.2,0.6,1)
task=c("medium")
#var.sigma =c(0.5,1,1.5)
var.sigma =c(3)

model<-"linear"
phi.beta<-NULL
intercept<-1
mean_y<-NULL

start=0


no<-2

n<-c(100,250,500)


for (i in delta0){
  for (ii in delta1){
    for (j in task){
      for (jj in no){
        if (jj==1) alteri=FALSE else alteri=TRUE

        for (z in n){
          for (zz in alteri){
            for (mm in var.sigma){
              nm<-paste("d0",i,"d1",ii,j,"a",ifelse(zz==TRUE,1,0),"no",jj,
                        "n",z,model,sep=".")

              if(model=="linear"|model=="quantile") nm<-paste(nm,"vs",mm,"in",intercept,sep=".") else
                nm<-paste(nm,"mY",mean_y,sep=".")

              if (model=="beta") nm<-paste(nm,"ph",phi.beta,sep=".")

              nm<-paste(nm,".txt",sep="")

              nmr1<-paste("results/res",nm,sep=".")

              dd<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(dd)!="try-error"){
                start=start+1



                my.sum<-function(x){
                  n<-length(x)
                  res<-c(n,
                         mean(x,na.rm=TRUE),
                         sd(x,na.rm=TRUE),
                         sum(is.na(x))/n*100
                  )
                  names(res)<-c("Nsim","mean","sd","pNA")
                  res

                }
                #res<-t(apply(dd[,-c(1:4)],2,my.sum))

                res.rej<-apply(dd[,1:4],2,function(x) mean(x<=0.05,na.rm=TRUE))
                res.nNA<-apply(dd[,1:4],2,function(x) sum(is.na(x)))

                #lmx<-100

                #cv1<-apply(dd1[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2<-apply(dd2[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv3<-apply(dd3[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)

                #cv1r<-apply(dd1r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2r<-apply(dd2r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)



                #make a dataframe with both methods in dif rows, add columns which will uniquely describe the setup
                dfi1<-as.data.frame(res.rej)
                dfi1$method<-c("GoF-S","GoF-P","F-S","F-P")

                dfi1$delta0<-i
                dfi1$delta1<-ii
                dfi1$task<-j
                dfi1$no<-jj
                dfi1$n<-z
                dfi1$alter<-zz
                dfi1$var.sigma<-mm
                dfi1$Nsim<-nrow(dd)
                dfi1$nNA<-res.nNA
                if (start==1) df<-dfi1 else df<-rbind(df,dfi1)
              }

            }}}}}
  }}

df_pvals<-df





methods_all<-c("F-P","F-S",
               "GoF-P","GoF-S")
methods_colors_all <- c(

  "F-P"       = "#009E73",  # green
  "F-S"        = "#006D5B",  # dark green

  "GoF-P"     = "#D55E00",  # vermillion
  "GoF-S"      = "#CC4452"   # red
)


method_labels_all<-c( "F-P","F-S",
                     "GoF-P","GoF-S")



which_display<-c(1:4)



methods_display<-methods_all[which_display]

method_colors <-methods_colors_all[which_display]


method_labels <- method_labels_all[which_display]


df_pvals$no<-2
df_pvals$method2<-df_pvals$method

df_pvals$h0<-2

df_pvals$rej05<-df_pvals$res.rej
######################

######################


#########layer 2
make_pval_plot <- function(methods = methods_display) {

  # ===========================================================================
  # DATA
  # ===========================================================================

  df_f <- df_pvals[
    df_pvals$no == 2 &
      df_pvals$h0 == df_pvals$no &
      df_pvals$method2 %in% methods,
  ]

  # Order sample sizes
  df_f$n <- factor(
    df_f$n,
    levels = c(100, 250, 500)
  )

  # Rejection rate as percentage
  df_f$rej05_pct <- 100 * df_f$rej05


  # ===========================================================================
  # SIMULATION MARGIN-OF-ERROR BAND
  # ===========================================================================

  p0 <- 0.05

  df_f$sim_moe <- 1.96 * sqrt(
    p0 * (1 - p0) / df_f$Nsim
  )

  df_f$band_low <- pmax(
    0,
    100 * (p0 - df_f$sim_moe)
  )

  df_f$band_high <- pmin(
    100,
    100 * (p0 + df_f$sim_moe)
  )

  # One gray band per sample-size facet
  band_df <- unique(
    df_f[, c("n", "band_low", "band_high")]
  )


  # ===========================================================================
  # MAIN PLOT
  # ===========================================================================

  p <- ggplot2::ggplot(
    df_f,
    ggplot2::aes(
      x = delta1,
      y = rej05_pct,
      color = method2,
      group = method2
    )
  ) +

    # Gray simulation margin-of-error band
    ggplot2::geom_rect(
      data = band_df,
      ggplot2::aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = band_low,
        ymax = band_high
      ),
      fill = "grey85",
      alpha = 1,
      color = NA,
      inherit.aes = FALSE
    ) +

    # Rejection-rate curves
    ggplot2::geom_line(
      alpha = 0.7
    ) +

    ggplot2::geom_point(
      alpha = 0.7,
      shape = 16
    ) +

    # Sample sizes in columns
    ggplot2::facet_grid(
      . ~ n,
      labeller = ggplot2::labeller(
        n = function(x) paste0("n = ", x)
      )
    ) +

    # Y axis
    #ggplot2::scale_y_continuous(
    #  breaks = c(0, 5, 25, 50, 75, 100),
    #  limits = c(0, 100)
    #) +

    # Method colors
    ggplot2::scale_color_manual(
      values = method_colors,
      labels = method_labels,
      name = "Method"
    ) +

    # Labels
    ggplot2::labs(
      x = expression(gamma[1]),
      y = "Rejection rate (%)"
    ) +

    # Theme
    ggplot2::theme_minimal(
      base_size = 14
    ) +

    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        color = "grey90"
      ),
      panel.grid.minor = ggplot2::element_blank(),

      strip.text = ggplot2::element_text(
        size = 13
      ),

      axis.title = ggplot2::element_text(
        size = 14
      ),

      axis.text = ggplot2::element_text(
        size = 11
      ),

      axis.text.x = ggplot2::element_text(
        size = 10
      ),

      legend.position = "right",
      legend.text = ggplot2::element_text(
        size = 12
      ),
      legend.title = ggplot2::element_text(
        size = 13
      )
    )

  return(p)
}


pdf("figs/rejection_setstart_2.pdf",height=8,width=8)
make_pval_plot()
dev.off()
#only show tab!



##table

library(dplyr)
library(tidyr)
library(xtable)

# Keep rows corresponding to the null hypothesis
df_ps <- df_pvals[df_pvals$no == df_pvals$h0, ]

# Calculate rejection rates, MCSEs, and percentages of missing p-values
df_ps <- df_ps %>%
  mutate(
    rej05_pct = 100 * rej05,
    mcse_pct = 100 * sqrt(rej05 * (1 - rej05) / Nsim),
    na_pct = 100 * nNA / Nsim,
    res = paste0(
      "\\parbox[t]{1.8cm}{\\centering ",
      sprintf("%.1f (%.1f)", rej05_pct, mcse_pct),
      "\\\\[-3pt]{\\tiny NA: ",
      sprintf("%.1f", na_pct),
      "}\\\\[3pt]}"
    )
  )

# Reshape results into a table
tab <- df_ps %>%
  select(h0, n, delta1, method2, res) %>%
  distinct() %>%
  pivot_wider(
    names_from = method2,
    values_from = res
  ) %>%
  arrange(h0, n, delta1)

# Column alignment (include the row-name column required by xtable)
align_vec <- c(
  "l",
  "r",
  "r",
  "r",
  rep("p{1.8cm}", ncol(tab) - 3)
)

# Generate LaTeX table
latex_tab <- print(
  xtable(
    tab,
    caption = paste(
      "Empirical rejection rates at $\\alpha=0.05$ (\\%),",
      "with Monte Carlo standard errors in parentheses and",
      "percentages of missing $p$-values for testing",
      "$H_0: K=K_0$ under the scenario with homoscedastic normal errors."
    ),
    label = "tab:results",
    align = align_vec,
    digits=c(0,0,0,1,2,2,2,2)
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity,
  print.results = FALSE
)



cat(latex_tab)


##main text table

library(dplyr)
library(tidyr)
library(xtable)

# Keep rows corresponding to the null hypothesis
df_ps <- df_pvals[df_pvals$no == df_pvals$h0, ]

# Calculate rejection rates and MCSEs
df_ps <- df_ps %>%
  mutate(
    rej05_pct = 100 * rej05,
    mcse_pct = 100 * sqrt(rej05 * (1 - rej05) / Nsim),
    res = ifelse(
      is.na(rej05_pct) | is.na(mcse_pct),
      "NA",
      sprintf("%.1f (%.1f)", rej05_pct, mcse_pct)
    )
  )

# Reshape results into a table
tab <- df_ps %>%
  select(h0, n, delta1, method2, res) %>%
  distinct() %>%
  pivot_wider(
    names_from = method2,
    values_from = res
  ) %>%
  arrange(h0, n, delta1)

# Column alignment
align_vec <- c(
  "l",
  "r",
  "r",
  "r",
  rep("c", ncol(tab) - 3)
)

# Generate LaTeX table
latex_tab_main <- print(
  xtable(
    tab,
    caption = paste(
      "Empirical rejection rates at $\\alpha=0.05$ (\\%),",
      "with Monte Carlo standard errors in parentheses for testing",
      "$H_0: K=K_0$ under the scenario with homoscedastic normal errors."
    ),
    label = "tab:results",
    align = align_vec,
    digits = c(0, 0, 0, 1, 2, 2, 2, 2)
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity,
  print.results = FALSE
)

cat(latex_tab_main)


plot(df_pvals$Nsim)




###############

##############some Fig that would explain this? Note we could also get p-values for the orig methods (same seeds): not the same no of iter though (yet)!
titles <- c(
  "Set position", "Set position",
  "Position not set", "Position not set"
)

nms1 <- c("All data", "GoF-S", "GoF-P", "F-S", "F-P")

rejection_rate <- function(pvals) {
  mean(pvals < 0.05, na.rm = TRUE) * 100
}

delta0 <- c(0.2)
delta1 <- c(0.2, 0.6, 1)
task <- c("medium")
var.sigma <- c(3)

model <- "linear"
phi.beta <- NULL
intercept <- 1
mean_y <- NULL

start <- 0
no <- 2
n <- c(100, 250, 500)

# Horizontal jitter width
jitter_width <- 0.15

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
                "a", ifelse(zz, 1, 0),
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

              nm <- paste0(nm, ".txt")

              nmr1 <- paste("results/res", nm, sep = ".")

              nmr2 <- paste(
                "D:/radulovic_joinpoint/simulation/numberOfChangepoints/revision_analysis/normal_hetero_residuals/results/pval",
                nm,
                sep = "."
              )

              dd <- try(
                read.table(
                  nmr1, sep = "\t",
                  fill = TRUE, header = FALSE
                ),
                silent = TRUE
              )

              if (!inherits(dd, "try-error")) {

                start <- start + 1

                dd1 <- try(
                  read.table(
                    nmr2, sep = "\t",
                    fill = TRUE, header = FALSE
                  ),
                  silent = TRUE
                )

                if (inherits(dd1, "try-error")) next

                names(dd1) <- paste0(
                  rep(c("gofsim", "gofperm", "fsim", "fperm"), each = 5),
                  0:4
                )

                # Relevant p-values when positions are not set
                dd1 <- dd1[, seq(from = 3, by = 5, length.out = 4)]

                # Match number of simulation replications
                nrw <- min(nrow(dd), nrow(dd1))

                dd <- dd[seq_len(nrw), , drop = FALSE]
                dd1 <- dd1[seq_len(nrw), , drop = FALSE]

                # Combine p-values
                df_pvs <- cbind(dd[, 1:4], dd1)

                names(df_pvs)[1:4] <- paste0(
                  c("gofsim", "gofperm", "fsim", "fperm"),
                  "setPosition"
                )

                # Estimated changepoint positions
                df_pos <- dd[, c(5, 6, 19, 20)]

                names(df_pos) <- c(
                  "tau1_setpos", "tau2_setpos",
                  "tau1", "tau2"
                )

                # Create output PDF
                pdf(
                  paste0("figs/fig", z, "_", ii, ".pdf"),
                  height = 8,
                  width = 10
                )

                # Plot order:
                # top row: position not set
                # bottom row: position set
                layout(matrix(c(3, 1, 4, 2), nrow = 2, byrow = TRUE))

                for (iii in 1:4) {

                  if (iii == 1 | iii == 3) {
                    ylabi <- expression(tau[1])
                    hlni <- 5
                  } else {
                    ylabi <- expression(tau[2])
                    hlni <- 10
                  }

                  # Relevant p-value columns
                  if (iii < 3) {
                    pcols <- 1:4
                  } else {
                    pcols <- 5:8
                  }

                  # Calculate rejection rates
                  rates <- sapply(
                    pcols,
                    function(k) rejection_rate(df_pvs[, k])
                  )

                  # Prepare boxplot data
                  box_data <- list(
                    df_pos[, iii],
                    df_pos[df_pvs[, pcols[1]] < 0.05, iii],
                    df_pos[df_pvs[, pcols[2]] < 0.05, iii],
                    df_pos[df_pvs[, pcols[3]] < 0.05, iii],
                    df_pos[df_pvs[, pcols[4]] < 0.05, iii]
                  )

                  # Initialize plotting region
                  plot(
                    NA,
                    xlim = c(0.5, 5.5),
                    ylim = c(0, 30),
                    xaxt = "n",
                    xlab = "",
                    ylab = ylabi,
                    main = titles[iii]
                  )

                  # Draw all observations FIRST (black points)
                  vals_all <- df_pos[, iii]
                  valid_all <- !is.na(vals_all)

                  points(
                    jitter(
                      rep(1, sum(valid_all)),
                      amount = jitter_width
                    ),
                    vals_all[valid_all],
                    pch = 16,
                    col = rgb(0, 0, 0, 0.4)
                  )

                  # Draw significant observations (red points)
                  for (kk in seq_along(pcols)) {

                    pcol <- pcols[kk]

                    selected <- !is.na(df_pvs[, pcol]) &
                      df_pvs[, pcol] < 0.05 &
                      !is.na(df_pos[, iii])

                    vals <- df_pos[selected, iii]

                    points(
                      jitter(
                        rep(kk + 1, length(vals)),
                        amount = jitter_width
                      ),
                      vals,
                      pch = 16,
                      col = rgb(1, 0, 0, 0.5)
                    )
                  }

                  # Draw boxplots ON TOP of points
                  boxplot(
                    box_data,
                    add = TRUE,
                    axes = FALSE,
                    #outline = TRUE,
                    boxwex = 0.5,
                    col = rgb(0.85, 0.85, 0.85, 0.7),
                    border = "black",
                    medcol = "black",
                    whiskcol = "black",
                    staplecol = "black",outline = FALSE
                  )

                  # Add x-axis labels
                  axis(
                    side = 1,
                    at = 1:5,
                    labels = nms1
                  )

                  # Add rejection rates below test labels
                  axis(
                    side = 1,
                    at = 2:5,
                    labels = paste0("(", round(rates, 1), "%)"),
                    line = 2,
                    tick = FALSE
                  )

                  # Add horizontal reference line last
                  abline(
                    h = hlni,
                    col = "gray40",
                    lty = 2,
                    lwd = 1
                  )
                }

                dev.off()
              }
            }
          }
        }
      }
    }
  }
}
