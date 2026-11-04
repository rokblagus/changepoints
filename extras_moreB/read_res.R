
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


no<-1

n<-c(100)


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

              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1
                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)


                names(ddd)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")
                names(dd1)<-c("gofsim1000","gofsim5000","gofsim10000")

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
                         sd((x - jj)^2,na.rm = TRUE)
                  )
                  names(res)<-c("Nsim","mean","sd","pNA","p0","p1","p2","p3","p4","p5","pCorrect","pUnder","pOver","sd_mse")
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

df_nchp=df



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


no<-1

n<-c(100)


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

              ddd<-try(read.table(nmdes,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)
              if (class(ddd)!="try-error"){
                start=start+1
                dd1<-try(read.table(nmr1,sep="\t",fill=TRUE,header=FALSE),silent=TRUE)


                names(ddd)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")
                names(dd1)<-paste0(rep(c("gofsim1000","gofsim5000","gofsim10000"),each=3),c(0:2))

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
                dfi1$method2<-rep(c("gofsim1000","gofsim5000","gofsim10000"),each=3)
                dfi1$h0<-rep(0:2,3)

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

df_pvals<-df


df_ps <- df_pvals %>%
  filter(no == h0) %>%
  mutate(
    rejection = paste0(
      round(rej05 * 100, 1),
      " (",
      round(sqrt(rej05 * (1 - rej05) / Nsim) * 100, 1),
      ")"
    )
  )


#------------------------------------------------------------
# Bias, MSE, % correct
#------------------------------------------------------------

df_nchp2 <- df_nchp %>%
  mutate(
    bias = mean - no,
    mcse_bias = sd / sqrt(Nsim),

    mse = bias^2 + sd^2,
    mcse_mse = sd_mse / sqrt(Nsim),

    mcse_pCorrect = sqrt(
      pCorrect * (100 - pCorrect) / Nsim
    ),

    bias_res = paste0(
      round(bias, 3),
      " (", round(mcse_bias, 3), ")"
    ),

    mse_res = paste0(
      round(mse, 3),
      " (", round(mcse_mse, 3), ")"
    ),

    correct_res = paste0(
      round(pCorrect, 1),
      " (", round(mcse_pCorrect, 1), ")"
    )
  ) %>%
  rename(method2 = method)


#------------------------------------------------------------
# Combine
#------------------------------------------------------------

tab_long <- df_ps %>%
  select(h0, n, delta1, method2, rejection) %>%
  left_join(
    df_nchp2 %>%
      select(
        no, n, delta1, method2,
        bias_res, mse_res, correct_res
      ),
    by = c(
      "h0" = "no",
      "n",
      "delta1",
      "method2"
    )
  ) %>%
  pivot_longer(
    cols = c(rejection, bias_res, mse_res, correct_res),
    names_to = "measure",
    values_to = "result"
  ) %>%
  mutate(
    measure = recode(
      measure,
      rejection = "Rejection 5\\%",
      bias_res = "Bias",
      mse_res = "MSE",
      correct_res = "\\% Correct"
    )
  )


tab <- tab_long %>%
  select(measure, delta1, method2, result) %>%
  pivot_wider(
    names_from = method2,
    values_from = result
  ) %>%
  arrange(
    delta1,
    factor(
      measure,
      levels = c(
        "Rejection 5\\%",
        "Bias",
        "MSE",
        "\\% Correct"
      )
    )
  )

print(
  xtable(
    tab,
    caption = "Rejection rates, bias, MSE, and percentage correct.",
    label = "tab:results",
    digits = 3
  ),
  include.rownames = FALSE,
  sanitize.text.function = identity
)

