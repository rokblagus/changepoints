

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


no<-1:3

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


                 nmr1<-paste("results/r1",nm,sep=".")

              dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(dd1)!="try-error"){
                start=start+1

                #make it conditional on the fact that the original model has converged - in LM this is irrelevant since lm always gives an estimate!


                my.sum<-function(x){
                  n<-length(x)
                  res<-c(n,
                         mean(x,na.rm=TRUE),
                         sd(x,na.rm=TRUE),
                         sum(is.na(x))/n*100,
                         sum(x==0,na.rm=TRUE)/n*100, #note that division is by nsim and not only those that converged!
                         sum(x==1,na.rm=TRUE)/n*100,
                         sum(x==2,na.rm=TRUE)/n*100,
                         sum(x==3,na.rm=TRUE)/n*100,
                         sum(x==4,na.rm=TRUE)/n*100,
                         sum(x==5,na.rm=TRUE)/n*100,
                         sum(x==jj,na.rm = TRUE)/n*100,
                         sum(x<jj,na.rm = TRUE)/n*100,
                         sum(x>jj,na.rm = TRUE)/n*100,
                         sd((x - jj)^2,na.rm = TRUE),
                         sum(x>(jj+1),na.rm = TRUE)/n*100
                  )
                  names(res)<-c("Nsim","mean","sd","pNA","p0","p1","p2","p3","p4","p5","pCorrect","pUnder","pOver","sd_mse","p_great1")
                  res

                }
                res<-t(apply(dd1[,1:2],2,my.sum))

                #lmx<-100

                #cv1<-apply(dd1[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2<-apply(dd2[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv3<-apply(dd3[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)

                #cv1r<-apply(dd1r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2r<-apply(dd2r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)



                #make a dataframe with both methods in dif rows, add columns which will uniquely describe the setup
                dfi1<-as.data.frame(res)
                dfi1$method<-c("Davies-F","Davies-R")

                dfi1$delta0<-i
                dfi1$delta1<-ii
                dfi1$task<-j
                dfi1$no<-jj
                dfi1$n<-z
                dfi1$alter<-zz
                dfi1$var.sigma<-mm

                if (start==1) df<-dfi1 else df<-rbind(df,dfi1)
              }

            }}}}}
  }}

dfalt<-df

###for null

delta0=c(0.2)
delta1=0
task=c("medium")
#var.sigma =c(0.5,1,1.5)
var.sigma =c(3)

model<-"linear"
phi.beta<-NULL
intercept<-1
mean_y<-NULL

start=0


no<-0

n<-c(100,250,500)


for (i in delta0){
  for (ii in delta1){
    for (j in task){
      for (jj in no){
        alteri=FALSE

        for (z in n){
          for (zz in alteri){
            for (mm in var.sigma){
              nm<-paste("d0",i,"d1",ii,j,"a",1,"no",jj,
                        "n",z,model,sep=".")

              if(model=="linear"|model=="quantile") nm<-paste(nm,"vs",mm,"in",intercept,sep=".") else
                nm<-paste(nm,"mY",mean_y,sep=".")

              if (model=="beta") nm<-paste(nm,"ph",phi.beta,sep=".")

              nm<-paste(nm,".txt",sep="")

                nmr1<-paste("results/r1",nm,sep=".")

              dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(dd1)!="try-error"){
                start=start+1



                my.sum<-function(x){
                  n<-length(x)
                  res<-c(n,
                         mean(x,na.rm=TRUE),
                         sd(x,na.rm=TRUE),
                         sum(is.na(x))/n*100,
                         sum(x==0,na.rm=TRUE)/n*100, #dont sum to 100, since n<n!
                         sum(x==1,na.rm=TRUE)/n*100,
                         sum(x==2,na.rm=TRUE)/n*100,
                         sum(x==3,na.rm=TRUE)/n*100,
                         sum(x==4,na.rm=TRUE)/n*100,
                         sum(x==5,na.rm=TRUE)/n*100,
                         sum(x==jj,na.rm = TRUE)/n*100,
                         sum(x<jj,na.rm = TRUE)/n*100,
                         sum(x>jj,na.rm = TRUE)/n*100,
                         sd((x - jj)^2,na.rm = TRUE),
                         sum(x>1,na.rm = TRUE)/n*100
                  )
                  names(res)<-c("Nsim","mean","sd","pNA","p0","p1","p2","p3","p4","p5","pCorrect","pUnder","pOver","sd_mse", "p_great1")
                  res

                }
                res<-t(apply(dd1[,1:2],2,my.sum))

                #lmx<-100

                #cv1<-apply(dd1[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2<-apply(dd2[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv3<-apply(dd3[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)

                #cv1r<-apply(dd1r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2r<-apply(dd2r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)



                #make a dataframe with both methods in dif rows, add columns which will uniquely describe the setup
                dfi1<-as.data.frame(res)
                dfi1$method<-c("Davies-F","Davies-R")

                dfi1$delta0<-i
                dfi1$delta1<-ii
                dfi1$task<-j
                dfi1$no<-jj
                dfi1$n<-z
                dfi1$alter<-zz
                dfi1$var.sigma<-mm

                if (start==1) df<-dfi1 else df<-rbind(df,dfi1)
              }

            }}}}}
  }}



library(dplyr)
library(tidyr)
library(ggplot2)

df_nchp_davies <- bind_rows(df, dfalt) #I guess this can be rbinded to all methods? yes, added in orig code


###plots only for these 2 methods

df_nchp <-df_nchp_davies

# Method labels
methods_all <- c(
  "Davies-F",
  "Davies-R"
)

method_labels_all <- c(
  "Davies-F",
  "Davies-R"
)

lookup <- setNames(method_labels_all, methods_all)

# Keep the original method variable for reference,
# but replace it with the plotting labels.
df_nchp <- df_nchp %>%
  mutate(
    method = recode(method, !!!lookup),
    method = factor(method, levels = method_labels_all)
  )

df_base <- df_nchp %>%
  mutate(
    pCorrect_prop = pCorrect / 100,
    pUnder_prop   = pUnder / 100,
    pOver_prop    = pOver / 100,
    pGreat1_prop  = p_great1 / 100
  )

# ---------------------------
# Bias and RMSE
# ---------------------------

df_bias_rmse <- bind_rows(

  # Bias = E(Khat) - K0
  df_base %>%
    mutate(
      measure = "Bias",
      estimate = mean - no,
      MCSE = sd / sqrt(Nsim)
    ),

  # RMSE = sqrt(Var(Khat) + Bias^2)
  df_base %>%
    mutate(
      measure = "RMSE",
      bias = mean - no,
      estimate = sqrt(sd^2 + bias^2),
      MCSE = sd_mse / (2 * sqrt(Nsim) * estimate)
    )
) %>%
  mutate(
    lower = estimate - 1.96 * MCSE,
    upper = estimate + 1.96 * MCSE,
    lower = if_else(measure == "RMSE",
                    pmax(0, lower), lower),
    measure = factor(measure, levels = c("Bias", "RMSE"))
  )

# ---------------------------
# Probability estimates
# ---------------------------

df_prob <- df_base %>%
  pivot_longer(
    cols = c(
      pCorrect_prop,
      pUnder_prop,
      pOver_prop,
      pGreat1_prop
    ),
    names_to = "measure",
    values_to = "p"
  ) %>%
  mutate(
    measure = recode(
      measure,
      pCorrect_prop = "Correct",
      pUnder_prop   = "Under",
      pOver_prop    = "Over",
      pGreat1_prop  = "Greater than 1"
    ),
    estimate = 100 * p,
    MCSE = 100 * sqrt(p * (1 - p) / Nsim),
    lower = pmax(0, estimate - 1.96 * MCSE),
    upper = pmin(100, estimate + 1.96 * MCSE)
  )

plot_results <- function(no_value, delta1_value = 0.6) {

  d_bias <- df_bias_rmse %>%
    filter(no == no_value)

  d_prob <- df_prob %>%
    filter(no == no_value)

  if (no_value > 0) {
    d_bias <- d_bias %>% filter(delta1 == delta1_value)
    d_prob <- d_prob %>% filter(delta1 == delta1_value)
  }
  # Select probabilities relevant to each case
  if (no_value == 0) {

    d_prob <- d_prob %>%
      filter(measure %in% c("Correct", "Over", "Greater than 1")) %>%
      mutate(
        measure = recode(
          measure,
          Correct = "P(hat(K) == 0)",
          Over = "P(hat(K) > 0)",
          `Greater than 1` = "P(hat(K) > 1)"
        )
      )

    prob_levels <- c(
      "P(hat(K) == 0)",
      "P(hat(K) > 0)",
      "P(hat(K) > 1)"
    )

  } else {

    d_prob <- d_prob %>%
      filter(measure %in% c("Correct", "Over", "Under")) %>%
      mutate(
        measure = recode(
          measure,
          Correct = "P(hat(K) == K[0])",
          Over = "P(hat(K) > K[0])",
          Under = "P(hat(K) < K[0])"
        )
      )

    prob_levels <- c(
      "P(hat(K) == K[0])",
      "P(hat(K) > K[0])",
      "P(hat(K) < K[0])"
    )
  }

  # Combine estimation-error and probability measures.
  # Bias and RMSE are on their original scale;
  # probabilities are expressed as percentages.
  d_bias <- d_bias %>%
    mutate(measure = as.character(measure))

  d_prob <- d_prob %>%
    mutate(measure = as.character(measure))

  d <- bind_rows(
    d_bias %>%
      select(method, n, measure, estimate, lower, upper),

    d_prob %>%
      select(method, n, measure, estimate, lower, upper)
  )

  measure_levels <- c("Bias", "RMSE", prob_levels)

  d <- d %>%
    mutate(
      measure = factor(measure, levels = measure_levels),
      n = factor(n, levels = sort(unique(n)))
    )

  ggplot(d, aes(x = method, y = estimate)) +

    # Reference line for zero bias
    geom_hline(
      data = d %>% filter(measure == "Bias"),
      aes(yintercept = 0),
      linetype = "dashed",
      color = "grey55",
      inherit.aes = FALSE
    ) +

    # 95% Monte Carlo confidence intervals
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.2,
      color = "grey50",
      linewidth = 0.5
    ) +

    geom_point(
      size = 2,
      color = "#2166AC"
    ) +

    # Metrics in rows, sample sizes in columns
    facet_grid(
      rows = vars(measure),
      cols = vars(n),
      scales = "free_y",
      labeller = label_parsed
    ) +

    labs(
      x = NULL,
      y = NULL
    )  +

    theme_bw() +

    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        size = 8
      ),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(size = 9),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}


p0 <- plot_results(0)
p1 <- plot_results(1,0.6)
p2 <- plot_results(2,0.6)
p3 <- plot_results(3,0.6)




p1 <- plot_results(1,0.2)
p2 <- plot_results(2,0.2)
p3 <- plot_results(3,0.2)



p1 <- plot_results(1,1)
p2 <- plot_results(2,1)
p3 <- plot_results(3,1)




####

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
  for (ii in c(1,2)){

    for (jj in no){


      for (z in n){

        for (mm in var.sigma){
          nm<-paste("d0",i,"type",ii,"a",0,"no",jj,
                    "n",z,model,sep=".")

          if(model=="linear"|model=="quantile") nm<-paste(nm,"vs",mm,"in",intercept,sep=".") else
            nm<-paste(nm,"mY",mean_y,sep=".")

          if (model=="beta") nm<-paste(nm,"ph",phi.beta,sep=".")

          nm<-paste(nm,".txt",sep="")


           nmr1<-paste("results_extras/r1",nm,sep=".")

           dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
           if (class(dd1)!="try-error"){
            start=start+1




            my.sum<-function(x){
              n<-length(x)
              res<-c(n,
                     mean(x,na.rm=TRUE),
                     sd(x,na.rm=TRUE),
                     sum(is.na(x))/n*100,
                     sum(x==0,na.rm=TRUE)/n*100, #note that division is by nsim and not only those that converged!
                     sum(x==1,na.rm=TRUE)/n*100,
                     sum(x==2,na.rm=TRUE)/n*100,
                     sum(x==3,na.rm=TRUE)/n*100,
                     sum(x==4,na.rm=TRUE)/n*100,
                     sum(x==5,na.rm=TRUE)/n*100,
                     sum(x==jj,na.rm = TRUE)/n*100,
                     sum(x<jj,na.rm = TRUE)/n*100,
                     sum(x>jj,na.rm = TRUE)/n*100,
                     sd((x - jj)^2,na.rm = TRUE),
                     sum(x>(jj+1),na.rm = TRUE)/n*100
              )
              names(res)<-c("Nsim","mean","sd","pNA","p0","p1","p2","p3","p4","p5","pCorrect","pUnder","pOver","sd_mse","p_great1")
              res

            }
            res<-t(apply(dd1[,1:2],2,my.sum))

            #lmx<-100

            #cv1<-apply(dd1[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
            #cv2<-apply(dd2[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
            #cv3<-apply(dd3[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)

            #cv1r<-apply(dd1r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
            #cv2r<-apply(dd2r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)



            #make a dataframe with both methods in dif rows, add columns which will uniquely describe the setup
            dfi1<-as.data.frame(res)
            dfi1$method<-c("Davies-F","Davies-R")

            dfi1$delta0<-i
            dfi1$type<-ii
            dfi1$no<-jj
            dfi1$n<-z
            dfi1$var.sigma<-mm

            if (start==1) df<-dfi1 else df<-rbind(df,dfi1)
          }

        }}}}}


dfalt<-df


df_nchp_davies <- bind_rows( dfalt)
