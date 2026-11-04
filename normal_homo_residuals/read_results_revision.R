 
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

              nmdes<-paste("results/d",nm,sep=".")
              nmr1<-paste("results/r1",nm,sep=".")
              if (z==500&jj==2&ii>0.5){
                nmdes2<-paste("results/mycompd",nm,sep=".")
                nmr12<-paste("results/mycompr1",nm,sep=".")

              }


              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1
                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                if (z==500&jj==2&ii>0.5){
                  ddd2<-try(read.table(nmdes2,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                  dd12<-try(read.table(nmr12,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)

                  if (nrow(ddd)>800) idr<-800 else idr<-nrow(ddd)

                  ddd<-rbind(ddd[1:idr,],ddd2)
                  dd1<-rbind(dd1[1:idr,],dd12)
                }

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
              nm<-paste("d0",i,"d1",ii,j,"a",ifelse(zz==TRUE,1,0),"no",jj,
                        "n",z,model,sep=".")

              if(model=="linear"|model=="quantile") nm<-paste(nm,"vs",mm,"in",intercept,sep=".") else
                nm<-paste(nm,"mY",mean_y,sep=".")

              if (model=="beta") nm<-paste(nm,"ph",phi.beta,sep=".")

              nm<-paste(nm,".txt",sep="")

              nmdes<-paste("results/d",nm,sep=".")
              nmr1<-paste("results/r1",nm,sep=".")

              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1
                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)


                names(ddd)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")
                names(dd1)<-c("BIC","AIC","BIC.sel","AIC.sel","davies","gofsim","gofperm","fsim","fperm",
                   "gofsimbonf","gofpermbonf","fsimbonf","fpermbonf")

                #cvf<-function(x,lim=1000) if(sum(ifelse(abs(x)>lim,1,0),na.rm=T)>0) 0 else 1

                #make it conditional on the fact that the original model has converged
                conv.or<-mean(ddd$convergence_init,na.rm=T)

                dd1<-dd1[ddd$convergence_init==1,]

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

 df_nchp <- bind_rows(df, dfalt)

 ##include Davies versions:

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


               nmr1<-paste("../different_davies/results/r1",nm,sep=".")

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

               nmr1<-paste("../different_davies/results/r1",nm,sep=".")

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




 df_nchp_davies <- bind_rows(df, dfalt) #I guess this can be rbinded to all methods?

 df_nchp<-rbind(df_nchp, df_nchp_davies)

 df_nchp<-df_nchp[-which(df_nchp$method=="davies"),]

 ###end include Davies versions


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


 ####add S methods that treat NAs differently

 #alt

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


               nmr1<-paste("../revised_S_methods/normal_hetero_residuals/results/r1",nm,sep=".")

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
               nm<-paste("d0",i,"d1",ii,j,"a",0,"no",jj,
                         "n",z,model,sep=".")

               if(model=="linear"|model=="quantile") nm<-paste(nm,"vs",mm,"in",intercept,sep=".") else
                 nm<-paste(nm,"mY",mean_y,sep=".")

               if (model=="beta") nm<-paste(nm,"ph",phi.beta,sep=".")

               nm<-paste(nm,".txt",sep="")

               nmr1<-paste("../revised_S_methods/normal_hetero_residuals/results/r1",nm,sep=".")

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

 df_nchp<-rbind(df_nchp,df,dfalt)

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


 p0 <- plot_results(0)
 p1 <- plot_results(1,0.6)
 p2 <- plot_results(2,0.6)
 p3 <- plot_results(3,0.6)

 pdf("figs/res_homo_0.pdf",height = 8,width=8)
 p0
 dev.off()

 pdf("figs/res_homo_1_06.pdf",height = 8,width=8)
 p1
 dev.off()
 pdf("figs/res_homo_2_06.pdf",height = 8,width=8)
 p2
 dev.off()
 pdf("figs/res_homo_3_06.pdf",height = 8,width=8)
 p3
 dev.off()


 p1 <- plot_results(1,0.2)
 p2 <- plot_results(2,0.2)
 p3 <- plot_results(3,0.2)

 pdf("figs/res_homo_1_02.pdf",height = 8,width=8)
 p1
 dev.off()
 pdf("figs/res_homo_2_02.pdf",height = 8,width=8)
 p2
 dev.off()
 pdf("figs/res_homo_3_02.pdf",height = 8,width=8)
 p3
 dev.off()


 p1 <- plot_results(1,1)
 p2 <- plot_results(2,1)
 p3 <- plot_results(3,1)

 pdf("figs/res_homo_1_1.pdf",height = 8,width=8)
 p1
 dev.off()
 pdf("figs/res_homo_2_1.pdf",height = 8,width=8)
 p2
 dev.off()
 pdf("figs/res_homo_3_1.pdf",height = 8,width=8)
 p3
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

              nmdes<-paste("results/d",nm,sep=".")
              nmr1<-paste("results/pval",nm,sep=".")
              if (z==500&jj==2&ii>0.5){
                nmdes2<-paste("results/mycompd",nm,sep=".")
                nmr12<-paste("results/mycomppval",nm,sep=".")


              }
              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1
                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                if (z==500&jj==2&ii>0.5){
                  ddd2<-try(read.table(nmdes2,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                  dd12<-try(read.table(nmr12,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)

                  if (nrow(ddd)>800) idr<-800 else idr<-nrow(ddd)

                  ddd<-rbind(ddd[1:idr,],ddd2)
                  dd1<-rbind(dd1[1:idr,],dd12)
                }

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

dfaltpv<-df


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
              nm<-paste("d0",i,"d1",ii,j,"a",ifelse(zz==TRUE,1,0),"no",jj,
                        "n",z,model,sep=".")

              if(model=="linear"|model=="quantile") nm<-paste(nm,"vs",mm,"in",intercept,sep=".") else
                nm<-paste(nm,"mY",mean_y,sep=".")

              if (model=="beta") nm<-paste(nm,"ph",phi.beta,sep=".")

              nm<-paste(nm,".txt",sep="")

              nmdes<-paste("results/d",nm,sep=".")
              nmr1<-paste("results/pval",nm,sep=".")

              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1
                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)


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


df_pvals<-rbind(df,dfaltpv)


##add S methods


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

              nmdes<-paste("..\\revised_S_methods\\normal_hetero_residuals/results/d",nm,sep=".")
              nmr1<-paste("..\\revised_S_methods\\normal_hetero_residuals/results/pval",nm,sep=".")

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

dfaltpv<-df


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
              nm<-paste("d0",i,"d1",ii,j,"a",ifelse(zz==TRUE,1,0),"no",jj,
                        "n",z,model,sep=".")

              if(model=="linear"|model=="quantile") nm<-paste(nm,"vs",mm,"in",intercept,sep=".") else
                nm<-paste(nm,"mY",mean_y,sep=".")

              if (model=="beta") nm<-paste(nm,"ph",phi.beta,sep=".")

              nm<-paste(nm,".txt",sep="")

              nmdes<-paste("..\\revised_S_methods\\normal_hetero_residuals/results/d",nm,sep=".")
              nmr1<-paste("..\\revised_S_methods\\normal_hetero_residuals/results/pval",nm,sep=".")

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

df_pvals<-rbind(df_pvals,df,dfaltpv)



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






######################

######################


#########layer 2

make_pval_plot <- function(no_val, methods = methods_display) {

  # ===========================================================================
  # DATA
  # ===========================================================================

  df_f <- df_pvals[
    df_pvals$no == no_val &
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

  if (no_val == 0) {

    # One band for each h0 column
    band_df <- unique(
      df_f[, c(
        "h0_label",
        "band_low",
        "band_high"
      )]
    )

  } else {

    # One band for each (n, h0) facet
    band_df <- unique(
      df_f[, c(
        "n",
        "h0_label",
        "band_low",
        "band_high"
      )]
    )
  }


  # ===========================================================================
  # BLUE RECTANGLE DATA
  # ===========================================================================

  if (no_val == 0) {

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

  } else {

    # Because we have n ~ h0_label facets, the rectangle
    # needs both facet variables.
    rect_df <- expand.grid(
      n = levels(df_f$n),
      h0_label = paste0(
        "H[0]: K == ",
        no_val
      ),
      stringsAsFactors = FALSE
    )

    rect_df$n <- factor(
      rect_df$n,
      levels = levels(df_f$n)
    )

    rect_df$h0_label <- factor(
      rect_df$h0_label,
      levels = paste0(
        "H[0]: K == ",
        h0_order
      )
    )
  }


  # ===========================================================================
  # PLOT: no = 0
  # ===========================================================================

  if (no_val == 0) {

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
    # PLOT: no >= 1
    # ===========================================================================

  } else {

    p <- ggplot2::ggplot(
      df_f,
      ggplot2::aes(
        x = delta1,
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

      # -----------------------------------------------------------------------
    # Columns = h0
    # Rows = n
    # -----------------------------------------------------------------------

    ggplot2::facet_grid(
      n ~ h0_label,
      labeller = ggplot2::labeller(
        h0_label = ggplot2::label_parsed,

        # IMPORTANT:
        # "=" rather than "==", so the displayed label is
        # n = 100, n = 250, etc.
        n = function(x) paste0(
          "n = ",
          x
        )
      )
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

        # Keep n labels horizontal
        strip.text.y = ggplot2::element_text(
          angle = 0
        ),

        legend.position = "none",

        axis.title = ggplot2::element_text(
          size = 14
        ),

        axis.text = ggplot2::element_text(
          size = 11
        ),

        axis.text.x = ggplot2::element_text(
          size = 10
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
      x = expression(gamma[1]),
      y = "Rejection rate (%)"
    )
  }


  # ===========================================================================
  # LEGEND
  # ===========================================================================

  p_leg <- ggplot2::ggplot(
    df_f,
    ggplot2::aes(
      x = delta1,
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


pdf("figs/rejection_homo_0.pdf",height=4,width=8)
make_pval_plot(0)
dev.off()


pdf("figs/rejection_homo_1.pdf",height=8,width=8)
make_pval_plot(1)
dev.off()


pdf("figs/rejection_homo_2.pdf",height=8,width=8)
make_pval_plot(2)
dev.off()



pdf("figs/rejection_homo_3.pdf",height=8,width=8)
make_pval_plot(3)
dev.off()


##table
library(dplyr)
library(tidyr)
library(xtable)

df_ps <- df_pvals %>%
  filter(no == h0) %>%
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

tab <- df_ps %>%
  select(h0, n, delta1, method2, res) %>%
  distinct() %>%
  pivot_wider(
    names_from = method2,
    values_from = res
  ) %>%
  arrange(h0, n, delta1)

# Column alignment (xtable requires an entry for the row-name column)
align_vec <- c(
  "l",
  "r",
  "r",
  "r",
  rep("p{1.8cm}", ncol(tab) - 3)
)

###here reorder!

tab<-tab[,c(1:3,5,4,8,7,6,9)]


# Generate LaTeX table code
latex_tab <- print(
  xtable(
    tab,
    caption = paste(
      "Empirical rejection rates, with Monte Carlo standard errors",
      "in parentheses and percentages of NA p-values on the second line."
    ),
    label = "tab:results",
    align = align_vec,
    digits = c(0, 0, 0, 1, 2, 2, 2, 2,2,2)
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity,
  print.results = FALSE
)


##for main paper/no NA

library(dplyr)
library(tidyr)
library(xtable)

df_ps <- df_pvals %>%
  filter(no == h0) %>%
  mutate(
    rej05_pct = 100 * rej05,
    mcse_pct = 100 * sqrt(rej05 * (1 - rej05) / Nsim),
    res = ifelse(
      is.na(rej05_pct) | is.na(mcse_pct),
      "NA",
      sprintf("%.1f (%.1f)", rej05_pct, mcse_pct)
    )
  )

tab <- df_ps %>%
  select(h0, n, delta1, method2, res) %>%
  distinct() %>%
  pivot_wider(
    names_from = method2,
    values_from = res
  ) %>%
  arrange(h0, n, delta1)

# Reorder columns
tab <- tab[, c(1:3, 5, 4, 8, 7, 6, 9)]

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
      "Empirical rejection rates, with Monte Carlo standard errors",
      "in parentheses."
    ),
    label = "tab:results",
    align = align_vec,
    digits = c(0, 0, 0, 1, 2, 2, 2, 2, 2, 2)
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity,
  print.results = FALSE
)



############################Layer 3

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

              nmdes<-paste("results/d",nm,sep=".")
              nmr1<-paste("results/times",nm,sep=".")
                nmr2<-paste("results/pval",nm,sep=".")

                if (z==500&jj==2&ii>0.5){
                  nmdes2<-paste("results/mycompd",nm,sep=".")
                  nmr12<-paste("results/mycomptimes",nm,sep=".")
                  nmr122<-paste("results/mycomppval",nm,sep=".")
                }

              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1
                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                 dd2<-try(read.table(nmr2,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                 if (z==500&jj==2&ii>0.5){
                   ddd2<-try(read.table(nmdes2,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                   dd12<-try(read.table(nmr12,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                   dd122<-try(read.table(nmr122,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                   if (nrow(ddd)>800) idr<-800 else idr<-nrow(ddd)

                   ddd<-rbind(ddd[1:idr,],ddd2)
                   dd1<-rbind(dd1[1:idr,],dd12)
                   dd2<-rbind(dd12[1:idr,],dd122)
                 }

                names(ddd)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")
                names(dd1)<-paste0(rep(c("gofsim","gofperm","fsim","fperm"),each=5),c(0:4))

                #cvf<-function(x,lim=1000) if(sum(ifelse(abs(x)>lim,1,0),na.rm=T)>0) 0 else 1

                #make it conditional on the fact that the original model has converged
                conv.or<-mean(ddd$convergence_init,na.rm=T)

                dd1<-dd1[ddd$convergence_init==1,]
                dd2<-dd2[ddd$convergence_init==1,]

                my.sum<-function(i,xx,yy){
x<-xx[,i]
y<-yy[,i]
                  n<-length(x)
                  na<-sum(is.na(x))
napv<-sum(is.na(y))
medianp<-median(x,na.rm=TRUE)
q3<-quantile(x,probs=0.75,na.rm=TRUE)
q1<-quantile(x,probs=0.25,na.rm=TRUE)
meanp<-mean(x,na.rm=TRUE)
sdp<-sd(x,na.rm=TRUE)

medianp2<-median(x[!is.na(y)],na.rm=TRUE)
q32<-quantile(x[!is.na(y)],probs=0.75,na.rm=TRUE)
q12<-quantile(x[!is.na(y)],probs=0.25,na.rm=TRUE)
meanp2<-mean(x[!is.na(y)],na.rm=TRUE)
sdp2<-sd(x[!is.na(y)],na.rm=TRUE)

res<-c(n,na,napv,meanp,sdp,medianp,q1,q3,meanp2,sdp2,medianp2,q12,q32)
                  names(res)<-c("Nsim","nNA","nNApv","meant","sdt","mediant","q1t","q3t","meant2","sdt2","mediant2","q1t2","q3t2")
                  res

                }
rr<-lapply(1:ncol(dd1),my.sum,dd1,dd2)
res<-matrix(unlist(rr),ncol=13,byrow=TRUE)
colnames(res)<-names(rr[[1]])
rownames(res)<-paste0(rep(c("gofsim","gofperm","fsim","fperm"),each=5),c(0:4))

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

dfalttimes<-df


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
              nm<-paste("d0",i,"d1",ii,j,"a",ifelse(zz==TRUE,1,0),"no",jj,
                        "n",z,model,sep=".")

              if(model=="linear"|model=="quantile") nm<-paste(nm,"vs",mm,"in",intercept,sep=".") else
                nm<-paste(nm,"mY",mean_y,sep=".")

              if (model=="beta") nm<-paste(nm,"ph",phi.beta,sep=".")

              nm<-paste(nm,".txt",sep="")

              nmdes<-paste("results/d",nm,sep=".")
              nmr1<-paste("results/times",nm,sep=".")
                       nmr2<-paste("results/pval",nm,sep=".")


              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1
                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)

                dd2<-try(read.table(nmr2,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
                names(ddd)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")
                 names(dd1)<-paste0(rep(c("gofsim","gofperm","fsim","fperm"),each=5),c(0:4))

                #cvf<-function(x,lim=1000) if(sum(ifelse(abs(x)>lim,1,0),na.rm=T)>0) 0 else 1

                #make it conditional on the fact that the original model has converged
                conv.or<-mean(ddd$convergence_init,na.rm=T)

                dd1<-dd1[ddd$convergence_init==1,]


   dd2<-dd2[ddd$convergence_init==1,]

                my.sum<-function(i,xx,yy){
x<-xx[,i]
y<-yy[,i]
                  n<-length(x)
                  na<-sum(is.na(x))
napv<-sum(is.na(y))
medianp<-median(x,na.rm=TRUE)
q3<-quantile(x,probs=0.75,na.rm=TRUE)
q1<-quantile(x,probs=0.25,na.rm=TRUE)
meanp<-mean(x,na.rm=TRUE)
sdp<-sd(x,na.rm=TRUE)

medianp2<-median(x[!is.na(y)],na.rm=TRUE)
q32<-quantile(x[!is.na(y)],probs=0.75,na.rm=TRUE)
q12<-quantile(x[!is.na(y)],probs=0.25,na.rm=TRUE)
meanp2<-mean(x[!is.na(y)],na.rm=TRUE)
sdp2<-sd(x[!is.na(y)],na.rm=TRUE)

res<-c(n,na,napv,meanp,sdp,medianp,q1,q3,meanp2,sdp2,medianp2,q12,q32)
                  names(res)<-c("Nsim","nNA","nNApv","meant","sdt","mediant","q1t","q3t","meant2","sdt2","mediant2","q1t2","q3t2")
                  res

                }
rr<-lapply(1:ncol(dd1),my.sum,dd1,dd2)
res<-matrix(unlist(rr),ncol=13,byrow=TRUE)
colnames(res)<-names(rr[[1]])
rownames(res)<-paste0(rep(c("gofsim","gofperm","fsim","fperm"),each=5),c(0:4))



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


df_times<-rbind(df,dfalttimes)

df_times_plot <- df_times %>%
  filter(
    h0 == no,
    no == 0 | delta1 == 0.6
  ) %>%
  mutate(
    method = lookup[method2],

    no = factor(
      no,
      levels = sort(unique(no))
    ),

    scenario = factor(
      paste0("n = ", n),
      levels = c(
        "n = 100",
        "n = 250",
        "n = 500"
      )
    ),

    nn = factor(
      paste0("K[0] == ", no),
      levels = c(
        "K[0] == 0",
        "K[0] == 1",
        "K[0] == 2",
        "K[0] == 3"
      )
    )
  )


plot_times <- ggplot(
  df_times_plot,
  aes(
    x = method,
    y = mediant2
  )
) +
  geom_errorbar(
    aes(
      ymin = q1t2,
      ymax = q3t2
    ),
    width = 0.2,
    color = "grey50",
    linewidth = 0.5
  ) +
  geom_point(
    size = 2,
    color = "#2166AC"
  ) +
  facet_grid(
    rows = vars(nn),
    cols = vars(scenario),
    labeller = labeller(
      nn = label_parsed
    )
  ) +
  scale_y_log10() +
  labs(
    x = NULL,
    y = "Median time (seconds)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),

    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 8
    ),

    strip.background = element_rect(
      fill = "grey95"
    ),

    strip.text = element_text(
      size = 8
    ),

    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

pdf("figs/times_homo.pdf",height=7,width=7)
plot_times
dev.off()









##make sense of everything: median accross scenarios by method

df_nchp$bias<-df_nchp$mean-df_nchp$no
df_nchp$rmse<-sqrt(df_nchp$bias**2+df_nchp$sd**2)

matrix(sapply(split(df_nchp$bias,list(df_nchp$method,df_nchp$no)),median),byrow=TRUE,nrow=4)

library(ggplot2)
library(dplyr)
library(tidyr)

# Exact order of methods
method_order <- c(
  "AIC", "AIC-S", "BIC", "BIC-S",
  "Davies-F", "Davies-R",
  "F-P", "F-P-B", "F-S", "F-S2", "F-S-B", "F-S2-B",
  "GoF-P", "GoF-P-B", "GoF-S", "GoF-S2", "GoF-S-B", "GoF-S2-B"
)

# Prepare data
df_base <- df_nchp %>%
  mutate(
    method = factor(method, levels = method_order),
    scenario = interaction(n, delta1, sep = ", delta1 = ")
  ) %>%
  filter(!is.na(method))

# Function to create a plot
make_plot <- function(data, metrics, metric_labels) {

  df_plot <- data %>%
    pivot_longer(
      cols = all_of(metrics),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = metrics,
        labels = unname(metric_labels[metrics])
      )
    )

  # Median across scenarios for each method, no, and metric
  median_df <- df_plot %>%
    group_by(metric, method, no) %>%
    summarise(
      median_value = median(value, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot(df_plot, aes(x = method, y = value, color = scenario)) +

    # Zero reference line for Bias only
    geom_hline(
      data = data.frame(
        metric = factor("Bias", levels = levels(df_plot$metric))
      ),
      aes(yintercept = 0),
      inherit.aes = FALSE,
      linetype = "dashed",
      color = "grey40"
    ) +

    # Individual scenario points
    geom_point(
      position = position_jitter(width = 0.1, height = 0),
      alpha = 1,
      size = 2
    ) +

    # Red squares for medians
    geom_point(
      data = median_df,
      aes(x = method, y = median_value),
      inherit.aes = FALSE,
      shape = 15,
      color = "red",
      size = 2
    ) +

    # Metrics in rows; K0 values in columns
    facet_grid(
      rows = vars(metric),
      cols = vars(no),
      scales = "free_y",
      labeller = labeller(
        metric = label_parsed,
        no = as_labeller(function(x) {
          paste0("K[0] == ", x)
        }, default = label_parsed)
      )
    ) +

    labs(
      x = "Method",
      y = NULL,
      color = expression("Scenario (" * n * ", " * delta[1] * ")")
    ) +

    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      strip.text.y = element_text(angle = 0),
      panel.grid.minor = element_blank()
    )
}

# ==================================================
# PLOT 1: K0 = 0
# Bias, RMSE, and P(K-hat = 0)
# ==================================================

metrics_zero <- c("bias", "rmse", "p0")

labels_zero <- c(
  bias = "Bias",
  rmse = "RMSE",
  p0 = "P(hat(K) == 0)"
)

df_zero <- df_base %>%
  filter(no == 0)

plot_K0_zero <- make_plot(
  data = df_zero,
  metrics = metrics_zero,
  metric_labels = labels_zero
)

print(plot_K0_zero)


# ==================================================
# PLOT 2: K0 > 0
# Bias, RMSE, P(K-hat = K0), P(K-hat > K0),
# and P(K-hat < K0)
# ==================================================

metrics_positive <- c(
  "bias", "rmse", "pCorrect", "pOver", "pUnder"
)

labels_positive <- c(
  bias = "Bias",
  rmse = "RMSE",
  pCorrect = "P(hat(K) == K[0])",
  pOver = "P(hat(K) > K[0])",
  pUnder = "P(hat(K) < K[0])"
)

df_positive <- df_base %>%
  filter(no > 0)

plot_K0_positive <- make_plot(
  data = df_positive,
  metrics = metrics_positive,
  metric_labels = labels_positive
)

print(plot_K0_positive)


###needs Figs, dataframe with nchp and description

#dir.create("results_for_analysis", showWarnings = FALSE)

# 1. Export the two figures as high-resolution PNG
ggsave(
  filename = "results_for_analysis/plot_K0_zero.png",
  plot = plot_K0_zero,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "results_for_analysis/plot_K0_positive.png",
  plot = plot_K0_positive,
  width = 18,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)

# 2. Export the underlying data
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


write.csv(
  df_times,
  file = "results_for_analysis/df_times.csv",
  row.names = FALSE
)



###convergence issues: (not yet completed!): two layers: estimating K and estimating p-values


#show as plot
df_pvals$pNA<-df_pvals$nNA/df_pvals$Nsim*100

#Layer 2 - p-values NAs


make_na_plot <- function(no_val, methods = methods_display) {

  # ===========================================================================
  # DATA
  # ===========================================================================

  df_f <- df_pvals[
    df_pvals$no == no_val &
      df_pvals$method2 %in% methods,
  ]

  # Order sample size
  df_f$n <- factor(
    df_f$n,
    levels = c(100, 250, 500)
  )

  # Percentage of simulation runs with NA p-values
  df_f$na_pct <- 100 * df_f$nNA / df_f$Nsim


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
  # BLUE RECTANGLE DATA: highlight the correct H0: K = no_val facet
  # ===========================================================================

  if (no_val == 0) {

    rect_df <- data.frame(
      h0_label = factor(
        paste0("H[0]: K == ", no_val),
        levels = paste0("H[0]: K == ", h0_order)
      )
    )

  } else {

    rect_df <- expand.grid(
      n = levels(df_f$n),
      h0_label = paste0("H[0]: K == ", no_val),
      stringsAsFactors = FALSE
    )

    rect_df$n <- factor(
      rect_df$n,
      levels = levels(df_f$n)
    )

    rect_df$h0_label <- factor(
      rect_df$h0_label,
      levels = paste0("H[0]: K == ", h0_order)
    )
  }


  # ===========================================================================
  # PLOT: no = 0
  # ===========================================================================

  if (no_val == 0) {

    p <- ggplot2::ggplot(
      df_f,
      ggplot2::aes(
        x = n,
        y = na_pct,
        color = method2,
        group = method2
      )
    ) +

      # Blue dashed rectangle around the H0: K = 0 facet
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

      # NA-proportion curves
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

      # Y axis
      ggplot2::scale_y_continuous(
        breaks = seq(0, 100, by = 10),
        limits = c(0, 100)
      ) +

      # Method colors
      ggplot2::scale_color_manual(
        values = method_colors,
        breaks = methods_display,
        labels = method_labels,
        name = "Method"
      ) +

      # Labels
      ggplot2::labs(
        x = "n",
        y = "Simulation runs with NA p-values (%)"
      )


    # ===========================================================================
    # PLOT: no >= 1
    # ===========================================================================

  } else {

    p <- ggplot2::ggplot(
      df_f,
      ggplot2::aes(
        x = delta1,
        y = na_pct,
        color = method2,
        group = method2
      )
    ) +

      # Blue dashed rectangle around the H0: K = no_val facet
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

      # NA-proportion curves
      ggplot2::geom_line(
        alpha = 0.7
      ) +

      ggplot2::geom_point(
        alpha = 0.7,
        shape = 16
      ) +

      # Columns = h0; rows = n
      ggplot2::facet_grid(
        n ~ h0_label,
        labeller = ggplot2::labeller(
          h0_label = ggplot2::label_parsed,
          n = function(x) paste0("n = ", x)
        )
      ) +

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

        strip.text.y = ggplot2::element_text(
          angle = 0
        ),

        legend.position = "none",

        axis.title = ggplot2::element_text(
          size = 14
        ),

        axis.text = ggplot2::element_text(
          size = 11
        ),

        axis.text.x = ggplot2::element_text(
          size = 10
        )
      ) +

      # Y axis
      ggplot2::scale_y_continuous(
        breaks = seq(0, 100, by = 10),
        limits = c(0, 100)
      ) +

      # Method colors
      ggplot2::scale_color_manual(
        values = method_colors,
        breaks = methods_display,
        labels = method_labels,
        name = "Method"
      ) +

      # Labels
      ggplot2::labs(
        x = expression(gamma[1]),
        y = "Simulation runs with NA p-values (%)"
      )
  }


  # ===========================================================================
  # LEGEND
  # ===========================================================================

  p_leg <- ggplot2::ggplot(
    df_f,
    ggplot2::aes(
      x = delta1,
      y = na_pct,
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


  legenda <- cowplot::get_legend(p_leg)


  # ===========================================================================
  # COMBINE PLOT + LEGEND
  # ===========================================================================

  cowplot::plot_grid(
    p,
    legenda,
    rel_widths = c(0.85, 0.15)
  )
}


# ===========================================================================
# CREATE AND SAVE PLOTS
# ===========================================================================

pdf("figs/NA_homo_0.pdf", height = 4, width = 8)
make_na_plot(0)
dev.off()

pdf("figs/NA_homo_1.pdf", height = 8, width = 8)
make_na_plot(1)
dev.off()

pdf("figs/NA_homo_2.pdf", height = 8, width = 8)
make_na_plot(2)
dev.off()

pdf("figs/NA_homo_3.pdf", height = 8, width = 8)
make_na_plot(3)
dev.off()


##check where we are
par(mfrow=c(1,3))
plot(df_nchp$Nsim)
plot(df_pvals$Nsim)
plot(df_times$Nsim)


###Table for H0: K=K0

cat(latex_tab) #supplementary table

#main file table (different, doesnt show nNA)

cat(latex_tab_main)


##final thing: NAs for estimating K

#pNA is percentage! Note: B methods have exactly the same rate of failures as non B, hence we could omit them?
#lets simply make a Table, it will not be too large
xt<-xtable(
  df_nchp[df_nchp$pNA>0,c(16,18,20,21,4)],
  digits=c(0,0,1,0,0,2)
)

print(xt,include.rownames=FALSE)
