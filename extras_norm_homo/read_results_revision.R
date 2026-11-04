

#this is for extra scenarios with different strength of changepoints, this is done!



###several levels of results: as before (layer 1), p-values (layer 2), times for refitted vs simulated F/GoF (layer 3)

####layer 1
 


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

              nmdes<-paste("results/d",nm,sep=".")
              nmr1<-paste("results/r1",nm,sep=".")


              nmdes2<-paste("results/mycompd",nm,sep=".")
              nmr12<-paste("results/mycompr1",nm,sep=".")

              nmdes12<-paste("results/2mycompd",nm,sep=".")
              nmr112<-paste("results/2mycompr1",nm,sep=".")


              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1

                ddd2<-try(read.table(nmdes2,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                ddd12<-try(read.table(nmdes12,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)

                if (nrow(ddd)<300) idr<-nrow(ddd) else idr<-300

                 ddd<-rbind(ddd[1:idr,],ddd12,ddd2)

                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                dd12<-try(read.table(nmr12,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                dd112<-try(read.table(nmr112,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)

                 dd1<-rbind(dd1[1:idr,],dd112,dd12)

                names(ddd)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")
                names(dd1)<-c("BIC","AIC","BIC.sel","AIC.sel","davies","gofsim","gofperm","fsim","fperm",
                   "gofsimbonf","gofpermbonf","fsimbonf","fpermbonf")

                #cvf<-function(x,lim=1000) if(sum(ifelse(abs(x)>lim,1,0),na.rm=T)>0) 0 else 1

                #make it conditional on the fact that the original model has converged - in LM this is irrelevant since lm always gives an estimate!
                conv.or<-mean(ddd$convergence_init,na.rm=T)

                dd1<-dd1[ddd$convergence_init==1,]

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
                  res<-t(apply(dd1,2,my.sum))

                #lmx<-100

                #cv1<-apply(dd1[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2<-apply(dd2[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv3<-apply(dd3[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)

                #cv1r<-apply(dd1r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2r<-apply(dd2r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)



                #make a dataframe with both methods in dif rows, add columns which will uniquely describe the setup
                dfi1<-as.data.frame(res)
                dfi1$method<-rownames(dfi1)

                dfi1$delta0<-i
                dfi1$type<-ii
                dfi1$no<-jj
                dfi1$n<-z
                dfi1$var.sigma<-mm

                if (start==1) df<-dfi1 else df<-rbind(df,dfi1)
              }

            }}}}}


dfalt<-df



 library(dplyr)
 library(tidyr)
 library(ggplot2)

 df_nchp <- bind_rows( dfalt)


 ###add Davies versions


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


           nmr1<-paste("../different_davies/results_extras/r1",nm,sep=".")

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


 df_nchp<-rbind(df_nchp, df_nchp_davies)

 df_nchp<-df_nchp[-which(df_nchp$method=="davies"),]




 ####end add Davies versions


 # Method labels
 methods_all <- c(
   "AIC", "AIC.sel", "BIC", "BIC.sel",  "Davies-F","Davies-R",
   "fperm", "fpermbonf", "fsim", "fsimbonf",
   "gofperm", "gofpermbonf", "gofsim", "gofsimbonf"
 )

 method_labels_all <- c(
   "AIC", "AIC-S", "BIC", "BIC-S",  "Davies-F","Davies-R",
   "F-P", "F-P-B", "F-S", "F-S-B",
   "GoF-P", "GoF-P-B", "GoF-S", "GoF-S-B"
 )

 lookup <- setNames(method_labels_all, methods_all)

 # Keep the original method variable for reference,
 # but replace it with the plotting labels.
 df_nchp <- df_nchp %>%
   mutate(
     method = recode(method, !!!lookup),
     method = factor(method, levels = method_labels_all)
   )



 ###ad S2 methods


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

           nmdes<-paste("..\\revised_S_methods\\extras_norm_homo/results/d",nm,sep=".")
           nmr1<-paste("..\\revised_S_methods\\extras_norm_homo/results/r1",nm,sep=".")





           ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
           if (class(ddd)!="try-error"){
             start=start+1


             dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)

             names(ddd)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")
             #names(dd1)<-c()

             #cvf<-function(x,lim=1000) if(sum(ifelse(abs(x)>lim,1,0),na.rm=T)>0) 0 else 1

             #make it conditional on the fact that the original model has converged - in LM this is irrelevant since lm always gives an estimate!
             conv.or<-mean(ddd$convergence_init,na.rm=T)

             dd1<-dd1[ddd$convergence_init==1,]

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
             res<-t(apply(dd1,2,my.sum))

             #lmx<-100

             #cv1<-apply(dd1[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
             #cv2<-apply(dd2[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
             #cv3<-apply(dd3[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)

             #cv1r<-apply(dd1r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
             #cv2r<-apply(dd2r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)



             #make a dataframe with both methods in dif rows, add columns which will uniquely describe the setup
             dfi1<-as.data.frame(res)
             dfi1$method<-c("GoF-S2","F-S2","GoF-S2-B","F-S2-B")

             dfi1$delta0<-i
             dfi1$type<-ii
             dfi1$no<-jj
             dfi1$n<-z
             dfi1$var.sigma<-mm

             if (start==1) df<-dfi1 else df<-rbind(df,dfi1)
           }

         }}}}}


 df_nchp<-rbind(df_nchp,df)

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

 plot_results <- function(no_value, type_value = 1) {

   d_bias <- df_bias_rmse %>%
     filter(no == no_value)

   d_prob <- df_prob %>%
     filter(no == no_value)

   if (no_value > 0) {
     d_bias <- d_bias %>% filter(type == type_value)
     d_prob <- d_prob %>% filter(type == type_value)
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

   method_order<-c("AIC","AIC-S","BIC","BIC-S","Davies-F","Davies-R","F-P","F-P-B","F-S","F-S2","F-S-B","F-S2-B",
                   "GoF-P", "GoF-P-B", "GoF-S",  "GoF-S2","GoF-S-B","GoF-S2-B")
   d <- d %>%
     mutate(
       method = factor(method, levels = method_order)
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


 p1 <- plot_results(2,1)
 p2 <- plot_results(2,2)



 pdf("figs/res_extras_1.pdf",height = 8,width=8)
 p1
 dev.off()
 pdf("figs/res_extras_2.pdf",height = 8,width=8)
 p2
 dev.off()



#######################Layer 2 - p-values


delta0=c(0.2)
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

              nmdes<-paste("results/d",nm,sep=".")
              nmr1<-paste("results/pval",nm,sep=".")


              nmdes2<-paste("results/mycompd",nm,sep=".")
              nmr12<-paste("results/mycomppval",nm,sep=".")

              nmdes12<-paste("results/2mycompd",nm,sep=".")
              nmr112<-paste("results/2mycomppval",nm,sep=".")



              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1


                ddd2<-try(read.table(nmdes2,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                ddd12<-try(read.table(nmdes12,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)

                if(nrow(ddd)<300) idr<-nrow(ddd) else idr<-300

                ddd<-rbind(ddd[1:idr,],ddd12,ddd2)

                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)

                dd12<-try(read.table(nmr12,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                dd112<-try(read.table(nmr112,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)

                dd1<-rbind(dd1[1:idr,],dd112,dd12)

                names(ddd)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")
                names(dd1)<-paste0(rep(c("gofsim","gofperm","fsim","fperm"),each=5),c(0:4))

                #cvf<-function(x,lim=1000) if(sum(ifelse(abs(x)>lim,1,0),na.rm=T)>0) 0 else 1

                #make it conditional on the fact that the original model has converged
                conv.or<-mean(ddd$convergence_init,na.rm=T)

                dd1<-dd1[ddd$convergence_init==1,]

                my.sum<-function(x){
                  n<-length(x)
                  na<-sum(is.na(x))
rej05<-mean(x<=0.05,na.rm=TRUE)
rej01<-mean(x<=0.01,na.rm=TRUE)
rej10<-mean(x<=0.1,na.rm=TRUE)
meanp<-mean(x,na.rm=TRUE)
sdp<-sd(x,na.rm=TRUE)
res<-c(n,na,rej05,rej01,rej10,meanp,sdp)
                  names(res)<-c("Nsim","nNA","rej05","rej01","rej10","meanp","sdp")
                  res

                }
                  res<-t(apply(dd1,2,my.sum))

                #lmx<-100

                #cv1<-apply(dd1[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2<-apply(dd2[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv3<-apply(dd3[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)

                #cv1r<-apply(dd1r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
                #cv2r<-apply(dd2r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)



                #make a dataframe with both methods in dif rows, add columns which will uniquely describe the setup
                dfi1<-as.data.frame(res)
                dfi1$method<-rownames(dfi1)
                  dfi1$method2<-rep(c("gofsim","gofperm","fsim","fperm"),each=5)
                dfi1$h0<-rep(0:4,4)

                dfi1$delta0<-i
                dfi1$type<-ii
                 dfi1$no<-jj
                dfi1$n<-z
                  dfi1$var.sigma<-mm

                if (start==1) df<-dfi1 else df<-rbind(df,dfi1)
              }

            }}}}}



dfaltpv<-df


df_pvals<-rbind(dfaltpv)

####add S2 methods

delta0=c(0.2)
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

          nmdes<-paste("..\\revised_S_methods\\extras_norm_homo/results/d",nm,sep=".")
          nmr1<-paste("..\\revised_S_methods\\extras_norm_homo/results/pval",nm,sep=".")






          ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
          if (class(ddd)!="try-error"){
            start=start+1



            dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)


            names(ddd)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")
            names(dd1)<-paste0(rep(c("GoF-S2","F-S2"),each=5),c(0:4))

            #cvf<-function(x,lim=1000) if(sum(ifelse(abs(x)>lim,1,0),na.rm=T)>0) 0 else 1

            #make it conditional on the fact that the original model has converged
            conv.or<-mean(ddd$convergence_init,na.rm=T)

            dd1<-dd1[ddd$convergence_init==1,]

            my.sum<-function(x){
              n<-length(x)
              na<-sum(is.na(x))
              rej05<-mean(x<=0.05,na.rm=TRUE)
              rej01<-mean(x<=0.01,na.rm=TRUE)
              rej10<-mean(x<=0.1,na.rm=TRUE)
              meanp<-mean(x,na.rm=TRUE)
              sdp<-sd(x,na.rm=TRUE)
              res<-c(n,na,rej05,rej01,rej10,meanp,sdp)
              names(res)<-c("Nsim","nNA","rej05","rej01","rej10","meanp","sdp")
              res

            }
            res<-t(apply(dd1,2,my.sum))

            #lmx<-100

            #cv1<-apply(dd1[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
            #cv2<-apply(dd2[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
            #cv3<-apply(dd3[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)

            #cv1r<-apply(dd1r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)
            #cv2r<-apply(dd2r[,c("x","U1.x","U2.x","U3.x","U4.x","U5.x")],1,cvf,lim=lmx)



            #make a dataframe with both methods in dif rows, add columns which will uniquely describe the setup
            dfi1<-as.data.frame(res)
            dfi1$method<-rownames(dfi1)
            dfi1$method2<-rep(c("GoF-S2","F-S2"),each=5)
            dfi1$h0<-rep(0:4,2)

            dfi1$delta0<-i
            dfi1$type<-ii
            dfi1$no<-jj
            dfi1$n<-z
            dfi1$var.sigma<-mm

            if (start==1) df<-dfi1 else df<-rbind(df,dfi1)
          }

        }}}}}


df_pvals<-rbind(df_pvals,df)


methods_all<-c("AIC","AIC.sel","BIC","BIC.sel","davies",
               "fperm","fpermbonf","fsim","fsimbonf",
               "gofperm","gofpermbonf","gofsim","gofsimbonf","F-S2","F-S2-B","GoF-S2","GoF-S2-B")
methods_colors_all <- c(
  "AIC"         = "#0072B2",  # blue
  "AIC.sel"     = "#56B4E9",  # light blue
  "BIC"         = "#E69F00",  # orange
  "BIC.sel"     = "#F0C75E",  # light orange
  "davies"      = "#CC79A7",  # purple/pink

  "fperm"       = "#009E73",  # green
  "fpermbonf"   = "#66C2A5",  # light green
  "fsim"        = "#006D5B",  # dark green
  "fsimbonf"    = "#A6D854",  # yellow-green


  "gofperm"     = "#D55E00",  # vermillion
  "gofpermbonf" = "#E69F7A", # light vermillion
  "gofsim"      = "#CC4452",  # red
  "gofsimbonf"  = "#F28E8E",  # light red

  "F-S2"        = "#007F8B",  # teal
  "F-S2-B"    = "#70C5C8",  # light teal

  "GoF-S2"      = "#8E5AA7",  # muted purple
  "GoF-S2-B"  = "#C5A0D5"
)


method_labels_all<-c("AIC","AIC-S","BIC","BIC-S","Davies",
                     "F-P","F-P-B","F-S","F-S-B",
                     "GoF-P","GoF-P-B","GoF-S","GoF-S-B","F-S2","F-S2-B","GoF-S2","GoF-S2-B")



which_display<-c(1:17)



methods_display<-methods_all[which_display]

method_colors <-methods_colors_all[which_display]


method_labels <- method_labels_all[which_display]






###convergence issues: only text

df_nchp[df_nchp$pNA>0,] #this is percentage!



######################

######################


#########layer 2

make_pval_plot <- function(type, methods = methods_display) {
no_val=2
  # ===========================================================================
  # DATA
  # ===========================================================================

  df_f <- df_pvals[
    df_pvals$type==type&
    df_pvals$no == 2 &
      df_pvals$method2 %in% methods,
  ]

  # Order n
  df_f$n <- factor(
    df_f$n,
    levels = c(100, 250, 500)
  )

  # Rejection rate as percentage
  df_f$rej05_pct <- 100 * df_f$rej05


  # ===========================================================================
  # FACET LABELS
  # ===========================================================================

  h0_order <- sort(unique(df_f$h0))

  df_f$h0_label <- paste0(
    "H[0]: K == ",
    df_f$h0
  )

  df_f$h0_label <- factor(
    df_f$h0_label,
    levels = paste0(
      "H[0]: K == ",
      h0_order
    )
  )


  # ===========================================================================
  # SIMULATION MARGIN OF ERROR
  #
  # Target rejection probability:
  #
  #     p = 0.05
  #
  # Margin of error:
  #
  #     1.96 * sqrt(0.05 * 0.95 / Nsim)
  #
  # Converted to percentages for the plot.
  # ===========================================================================

  p0 <- 0.05

  df_f$sim_moe <- 1.96 * sqrt(
    p0 * (1 - p0) / df_f$Nsim
  )

  df_f$band_low <- 100 * (
    p0 - df_f$sim_moe
  )

  df_f$band_high <- 100 * (
    p0 + df_f$sim_moe
  )

  # Keep the band within the plotting range
  df_f$band_low <- pmax(
    df_f$band_low,
    0
  )

  df_f$band_high <- pmin(
    df_f$band_high,
    100
  )


  # ===========================================================================
  # GRAY BAND DATA
  # ===========================================================================



    # One band for each h0 column
    band_df <- unique(
      df_f[, c(
        "h0_label",
        "band_low",
        "band_high"
      )]
    )



  # ===========================================================================
  # BLUE RECTANGLE DATA
  # ===========================================================================


    # One rectangle around the correct h0 column
    rect_df <- data.frame(
      h0_label = factor(
        paste0(
          "H[0]: K == ",
          no_val
        ),
        levels = paste0(
          "H[0]: K == ",
          h0_order
        )
      )
    )


  # ===========================================================================
  # PLOT: no = 0
  # ===========================================================================

  p <- ggplot2::ggplot(
    df_f,
    ggplot2::aes(
      x = n,
      y = rej05_pct,
      color = method2,
      group = method2
    )
  ) +

    # -----------------------------------------------------------------------
  # Uniform light-gray simulation band
  # -----------------------------------------------------------------------

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

    # -----------------------------------------------------------------------
  # Blue dashed rectangle around correct H0:K=no_val column
  # -----------------------------------------------------------------------

  ggplot2::geom_rect(
    data = rect_df,
    ggplot2::aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "blue",
    alpha = 0.05,
    color = "blue",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +

    # -----------------------------------------------------------------------
  # Rejection-rate curves
  # -----------------------------------------------------------------------

  ggplot2::geom_line(
    alpha = 0.7
  ) +

    ggplot2::geom_point(
      alpha = 0.7,
      shape = 16
    ) +

    # Columns = h0
    ggplot2::facet_grid(
      . ~ h0_label,
      labeller = ggplot2::label_parsed
    ) +

    # -----------------------------------------------------------------------
  # Theme
  # -----------------------------------------------------------------------

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

      legend.position = "none",

      axis.title = ggplot2::element_text(
        size = 14
      ),

      axis.text = ggplot2::element_text(
        size = 11
      ),

      axis.text.x = ggplot2::element_text(
        size = 10,
        angle = 45,
        hjust = 1
      )
    ) +

    # -----------------------------------------------------------------------
  # Y axis
  # -----------------------------------------------------------------------

  ggplot2::scale_y_continuous(
    breaks = c(
      0,
      5,
      25,
      50,
      75,
      100
    ),
    limits = c(
      0,
      100
    )
  ) +

    # -----------------------------------------------------------------------
  # Method colors
  # -----------------------------------------------------------------------

  ggplot2::scale_color_manual(
    values = method_colors,
    breaks = methods_display,
    labels = method_labels,
    name = "Method"
  ) +

    # -----------------------------------------------------------------------
  # Labels
  # -----------------------------------------------------------------------

  ggplot2::labs(
    x = "n",
    y = "Rejection rate (%)"
  )



  # ===========================================================================
  # LEGEND
  # ===========================================================================

  p_leg <- ggplot2::ggplot(
    df_f,
    ggplot2::aes(
      x = type,
      y = rej05_pct,
      color = method2,
      group = method2
    )
  ) +

    ggplot2::geom_line() +

    ggplot2::geom_point() +

    ggplot2::scale_color_manual(
      values = method_colors,
      breaks = methods_display,
      labels = method_labels,
      name = "Method"
    ) +

    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(
        size = 12
      ),
      legend.title = ggplot2::element_text(
        size = 13
      )
    )


  legenda <- cowplot::get_legend(
    p_leg
  )


  # ===========================================================================
  # COMBINE PLOT + LEGEND
  # ===========================================================================

  cowplot::plot_grid(
    p,
    legenda,
    rel_widths = c(
      0.85,
      0.15
    )
  )
}



methods_display <- c(
  "fperm", "fsim", "F-S2",
  "gofperm", "gofsim", "GoF-S2",
  "AIC", "AIC.sel", "BIC", "BIC.sel", "davies",
  "fpermbonf", "fsimbonf", "gofpermbonf", "gofsimbonf",
  "F-S2-B", "GoF-S2-B"
)

# Colors matched by method name
method_colors <- methods_colors_all[methods_display]

# Create a named lookup: method name -> correct display label
labels_lookup <- setNames(method_labels_all, methods_all)

# Retrieve labels in exactly the same order as methods_display
method_labels <- labels_lookup[methods_display]

# Keep names aligned with the method names for ggplot2
names(method_labels) <- methods_display

pdf("figs/rejection_extra_1.pdf",height=4,width=8)
make_pval_plot(1)
dev.off()


pdf("figs/rejection_extra_2.pdf",height=4,width=8)
make_pval_plot(2)
dev.off()



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
  select(h0, n, type, method2, res) %>%
  distinct() %>%
  pivot_wider(
    names_from = method2,
    values_from = res
  ) %>%
  arrange(h0, n, type)

# Remove h0 column, as in the original code
tab <- tab[, -1]

# Recode type
tab$type <- ifelse(tab$type == 1, 10, 25)

# Preserve the original column order
tab <- tab[, c(1:2, 4, 3, 7, 6, 5, 8)]

# Set column alignment (include the row-name column required by xtable)
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
    digits=c( 0, 0, 0, 2, 2, 2, 2,2,2)
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity,
  print.results = FALSE
)



cat(latex_tab)




df_ps<-df_pvals[df_pvals$no==df_pvals$h0,]

df_ps$res<-paste0(round(df_ps$rej05*100,1), " (",round(sqrt(df_ps$rej05*(1-df_ps$rej05)/df_ps$Nsim)*100,1),")")

library(dplyr)
library(tidyr)
library(xtable)

tab <- df_ps %>%
  mutate(
    scenario = paste0(
      "K_0=", h0,
      ", n=", n,
      ", type=", type
    )
  ) %>%
  select(scenario, method2, res) %>%
  distinct() %>%
  pivot_wider(
    names_from = method2,
    values_from = res
  )

tab <- df_ps %>%
  select(h0, n, type, method2, res) %>%
  distinct() %>%
  pivot_wider(
    names_from = method2,
    values_from = res
  ) %>%
  arrange(h0, n, type)

tab<-tab[,-1]
tab$type<-ifelse(tab$type==1,10,25)

tab<-tab[,c(1:2,4,3,7,6,5,8)]


library(xtable)
print(
  xtable(
    tab,
    caption = "Results by null hypothesis, sample size, and effect size.",
    label = "tab:results",
    digits = c( 0, 0, 0, 2, 2, 2, 2,2,2)
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity
)

par(mfrow=c(1,2))
plot(df_nchp$Nsim)
plot(df_pvals$Nsim)



df_nchp$bias<-df_nchp$mean-df_nchp$no
df_nchp$rmse<-sqrt(df_nchp$bias**2+df_nchp$sd**2)

dir.create("results_for_analysis", showWarnings = FALSE)



write.csv(
  df_nchp,
  file = "results_for_analysis/df_nchp.csv",
  row.names = FALSE
)


write.csv(
  df_pvals,
  file = "results_for_analysis/df_pvals.csv",
  row.names = FALSE
)



#####convergence issues


#layer 2 p-values


# ============================================================
# Plot percentage of NA estimates
# Uses df_pvals created by the supplied simulation code
# ============================================================

library(ggplot2)
library(cowplot)

# Methods present in the supplied code
methods_na <- c(
  "fperm", "fsim", "F-S2",
  "gofperm", "gofsim", "GoF-S2"
)

# Match the colors and display labels defined in your script
method_colors_na <- methods_colors_all[methods_na]

labels_lookup <- setNames(method_labels_all, methods_all)
method_labels_na <- labels_lookup[methods_na]
names(method_labels_na) <- methods_na


make_na_plot <- function(type, methods = methods_na) {

  no_val <- 2

  # ----------------------------------------------------------
  # DATA
  # ----------------------------------------------------------

  df_f <- df_pvals[
    df_pvals$type == type &
      df_pvals$no == no_val &
      df_pvals$method2 %in% methods,
  ]

  if (nrow(df_f) == 0) {
    stop("No data found for this type and method selection.")
  }

  # Order sample sizes
  df_f$n <- factor(
    df_f$n,
    levels = c(100, 250, 500)
  )

  # Percentage of simulation runs with NA estimates
  df_f$na_pct <- 100 * df_f$nNA / df_f$Nsim

  # ----------------------------------------------------------
  # FACET LABELS
  # ----------------------------------------------------------

  h0_order <- sort(unique(df_f$h0))

  df_f$h0_label <- factor(
    paste0("H[0]: K == ", df_f$h0),
    levels = paste0("H[0]: K == ", h0_order)
  )

  # ----------------------------------------------------------
  # BLUE RECTANGLE AROUND H0: K = no_val
  # ----------------------------------------------------------

  rect_df <- data.frame(
    h0_label = factor(
      paste0("H[0]: K == ", no_val),
      levels = paste0("H[0]: K == ", h0_order)
    )
  )

  # ----------------------------------------------------------
  # MAIN PLOT
  # ----------------------------------------------------------

  p <- ggplot2::ggplot(
    df_f,
    ggplot2::aes(
      x = n,
      y = na_pct,
      color = method2,
      group = method2
    )
  ) +

    # Highlight the correctly specified H0: K = 2 panel
    ggplot2::geom_rect(
      data = rect_df,
      ggplot2::aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "blue",
      alpha = 0.05,
      color = "blue",
      linetype = "dashed",
      inherit.aes = FALSE
    ) +

    # NA-percentage curves
    ggplot2::geom_line(
      alpha = 0.7
    ) +

    ggplot2::geom_point(
      alpha = 0.7,
      shape = 16
    ) +

    # One panel for each value of H0: K
    ggplot2::facet_grid(
      . ~ h0_label,
      labeller = ggplot2::label_parsed
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
      legend.position = "none",
      axis.title = ggplot2::element_text(
        size = 14
      ),
      axis.text = ggplot2::element_text(
        size = 11
      ),
      axis.text.x = ggplot2::element_text(
        size = 10,
        angle = 45,
        hjust = 1
      )
    ) +

    # Y-axis in percentages
    ggplot2::scale_y_continuous(
      breaks = seq(0, 100, by = 10),
      limits = c(0, 100)
    ) +

    # Method colors and labels
    ggplot2::scale_color_manual(
      values = method_colors_na,
      breaks = methods,
      labels = method_labels_na[methods],
      name = "Method"
    ) +

    # Axis labels
    ggplot2::labs(
      x = "n",
      y = "Simulation runs with NA estimates (%)"
    )

  # ----------------------------------------------------------
  # LEGEND
  # ----------------------------------------------------------

  p_leg <- ggplot2::ggplot(
    df_f,
    ggplot2::aes(
      x = n,
      y = na_pct,
      color = method2,
      group = method2
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(
      values = method_colors_na,
      breaks = methods,
      labels = method_labels_na[methods],
      name = "Method"
    ) +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 13)
    )

  legenda <- cowplot::get_legend(p_leg)

  # Combine plot and legend
  cowplot::plot_grid(
    p,
    legenda,
    rel_widths = c(0.85, 0.15)
  )
}


# ============================================================
# SAVE PLOTS FOR BOTH SIMULATION TYPES
# ============================================================


pdf("figs/NA_estimates_type1.pdf", height = 4, width = 10)
make_na_plot(type = 1)
dev.off()

pdf("figs/NA_estimates_type2.pdf", height = 4, width = 10)
make_na_plot(type = 2)
dev.off()


##final thing: NAs for estimating K

#pNA is percentage! Note: B methods have exactly the same rate of failures as non B, hence we could omit them?
#lets simply make a Table, it will not be too large

  df_nchp[df_nchp$pNA>0 ,] #no issues in this example ?

