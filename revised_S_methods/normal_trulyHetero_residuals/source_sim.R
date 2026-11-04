

###required packages


library(mvtnorm)
library(purrr)
library(brglm2)
library(segmented)
library(data.table)

setDTthreads(1)
###segmented does not work correctly with bias reduced estimator!, it works but you have to use segmented.default!

###source my segmented hack (to force him to use BR estimators!), not needed (thank god!)

#source("rok_segmented_glm.R")


###needed so that segmented.default works also for QR (there are different options for geting se, I am using ker, see vignete for details)

vcov.rq<-function(x){
  sm<-summary(x,se = "ker")$coef[,2]
  mat<-diag(sm**2)
  colnames(mat)<-rownames(mat)<-names(sm)
  mat
}

CvM<-function(x){ sum( x**2  )  }

########sim Z data

#the same as in Geroldinger et al. 2022 (see paper and supplement for reference)
#from paper:

#For the generation of the five covariates X1, ?X5 we followed ideas by Binder et al. [22] in order to obtain realistic data. By first sampling from five standard normally distributed variables Z1, ?Z5 with correlation structure as defined in Table S2 and then applying the transformations specified in Table S2 we obtained three binary variables X1, X2, X3, one ordinal variable X4 and one continuous variable X5. The covariates X1 and X2 were generated as ?between-cluster? covariates by requiring Z1 and Z2 to be constant within clusters, while the three other covariates were generated as ?within-cluster? covariates, i.e. they were allowed to vary between the observations of a cluster. To avoid extreme values in the metric covariate X5, we winsorized it at the third quartile plus three times the interquartile range and at the first quartile minus three times the interquartile range, where the quartiles were determined by applying the transformation given in Table S2 to the quartiles of the standard normally distributed variable Z5.
#Reference: Binder H, Sauerbrei W, Royston P. Multivariable model-building with continuous covariates: 1. Performance measures and simulation design. Technical report FDM-preprint 105. 2011; University of Freiburg, Germany.

simZ<-function(n){

  cor.mat<-rbind(c(1,0.5,0.5,0,0),
                 c(0.5,1,0,0.5,0),
                 c(0.5,0,1,0,-0.3),
                 c(0,0.5,0,1,0.5),
                 c(0,0,-0.3,0.5,1))
  X<-rmvnorm(n,sigma=cor.mat,method="chol")
  Z1<-ifelse(X[,1]<(-1),1,0)
  Z2<-ifelse(X[,2]<(-1),1,0)
  Z3<-ifelse(X[,3]<0,1,0)
  Z4<-ifelse(X[,4]>=0.5,1,0)+ifelse(X[,4]>=1.5,1,0)
  Z5<-10*X[,5]+55
  upl<-quantile(X[,5],probs=0.75)+3*(quantile(X[,5],probs=0.75)-quantile(X[,5],probs=0.25))
  lll<-quantile(X[,5],probs=0.75)-3*(quantile(X[,5],probs=0.75)-quantile(X[,5],probs=0.25))

  upll<-upl*10+55
  llll<-lll*10+55
  Z5[Z5<llll]<-llll
  Z5[Z5>upll]<-upll

  cbind(Z1,Z2,Z3,Z4,Z5)
  #Z<-model.matrix(~Z1+Z2+Z3+factor(Z4)+Z5-1)
  #Z
}


##other aux fncts to sim data

#pos of chpnts
get.tau<-function(a=0,b=30,task,no){

  if (task=="easy") f<-1

  if (task=="medium") f<-2

  if (task=="hard") f<-4

  floor(seq(from=a,to=b/f,length.out = 2+no)[-c(1,2+no)])

}


#sim all data needed to fit the model

#argument var.sigma only applies to linear and quantile
#argument phi.beta only applies to beta, note var(Y)=mu(1-mu)/(1+phi)

mysim_data<-function(n,dif,no, delta0,delta1,alter,model,intercept,var.sigma,phi.beta){

  x<-rdunif(n, b=30 , a=0)# x must be greater than 30


  Z<-simZ(n)
  betas<-c(0.69,0.69,0.69,0.35,-0.035)
  lp1<-Z%*%betas


  if(no==0){
    lp2<-intercept+delta0*x
    #lp<-lp2 +lp1
    chpoints=NA
    # endif
    # return("easyvvvvv")
  }
  else {
    tau<-get.tau(a=0,b=30,dif,no)
    chpoints=tau
    if(no==1)
    {

      tau1<-tau[1]
      lp2<-intercept+delta0*x+delta1* ifelse(x-tau1<0,0, x-tau1) #1 chpnt +epsilon
      #lp<-lp2 +lp1




    }

    if(no==2)
    {
      tau1=tau[1]
      tau2=tau[2]
      if (alter==FALSE) lp2<-intercept+delta0*x+delta1* ifelse(x-tau1<0,0, x-tau1)+delta1* ifelse(x-tau2<0,0, x-tau2)
      else
        lp2<-intercept+delta0*x+delta1* ifelse(x-tau1<0,0, x-tau1)-delta1* ifelse(x-tau2<0,0, x-tau2)





    }


    if(no==3)
    {
      tau1=tau[1]
      tau2=tau[2]
      tau3=tau[3]

      if (alter==FALSE) lp2<-intercept+delta0*x+delta1* ifelse(x-tau1<0,0, x-tau1)+delta1* ifelse(x-tau2<0,0, x-tau2)+delta1* ifelse(x-tau3<0,0, x-tau3) #1 chpnt
      else
        lp2<-intercept+delta0*x+delta1* ifelse(x-tau1<0,0, x-tau1)-delta1* ifelse(x-tau2<0,0, x-tau2)+delta1* ifelse(x-tau3<0,0, x-tau3)


    }


    if(no==4)
    {
      tau1=tau[1]
      tau2=tau[2]
      tau3=tau[3]
      tau4=tau[4]
      if (alter==FALSE) lp2<-intercept+delta0*x+delta1* ifelse(x-tau1<0,0, x-tau1)+delta1* ifelse(x-tau2<0,0, x-tau2)+delta1* ifelse(x-tau3<0,0, x-tau3)+
        delta1* ifelse(x-tau4<0,0, x-tau4)
      else
        lp2<-intercept+delta0*x+delta1* ifelse(x-tau1<0,0, x-tau1)-delta1* ifelse(x-tau2<0,0, x-tau2)+delta1* ifelse(x-tau3<0,0, x-tau3)-
        delta1* ifelse(x-tau4<0,0, x-tau4)
    }


    if(no==5)
    {
      tau1=tau[1]
      tau2=tau[2]
      tau3=tau[3]
      tau4=tau[4]
      tau5=tau[5]
      if (alter==FALSE) lp2<-intercept+delta0*x+delta1* ifelse(x-tau1<0,0, x-tau1)+delta1* ifelse(x-tau2<0,0, x-tau2)+delta1* ifelse(x-tau3<0,0, x-tau3)+
        delta1* ifelse(x-tau4<0,0, x-tau4)+    delta1* ifelse(x-tau5<0,0, x-tau5)
      else
        lp2<-intercept+delta0*x+delta1* ifelse(x-tau1<0,0, x-tau1)-delta1* ifelse(x-tau2<0,0, x-tau2)+delta1* ifelse(x-tau3<0,0, x-tau3)-
        delta1* ifelse(x-tau4<0,0, x-tau4)+    delta1* ifelse(x-tau5<0,0, x-tau5)

    }









    #return("easy")





  }


  lp<-lp2 +lp1

  if (model=="linear"|model=="quantile"){
       var.sigma.het<-(1+0.25*x)
  eps<-rnorm(n,sd=sqrt(var.sigma.het))
    y<-lp +eps
  }
  if (model=="logistic"){
    prob<-1/(1+exp(-lp))
    y<-rbinom(n,prob=prob,size=1)
  }
  if (model=="poisson"){
    rate<-exp(lp)

    y<-rpois(n,lambda=rate)

  }
  if (model=="beta"){
    mu<-1/(1+exp(-lp))
    #phi<-phi.beta #should not be too small!
    var<-mu*(1-mu)/(1+phi.beta)
    estBetaParams <- function(mu, var) {
      alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
      beta <- alpha * (1 / mu - 1)
      return(params = list(alpha = alpha, beta = beta))
    }
    sp<-estBetaParams(mu,var)

    #shape1<-mu*phi
    #shape2<-phi-shape1
    shape1<-sp$alpha
    shape2<-sp$beta
    y<-rbeta(n,shape1,shape2)
    y[y<1e-8]<-1e-8
    y[y>1-1e-8]<-1-1e-8
  }

  df<-data.frame(Y=y,Z =Z , X = x)

  my_list <- list(data=df , "chpoints"=chpoints)
  return(my_list)



}



#fntcs related to finding the intercept for different models (not relevant for linear and quantile)

#poisson

get.beta0_pois<-function(beta0, lp, target.mean){
  -target.mean + mean(exp(beta0 + lp))
}

get_intercept_poisson<-function(mn,delta0,delta1,chpnts,alter){
  #lpZ<-c(0.69,0.69,0.69,0.35,-0.035)%*%matrix(c(0.158,0.158,0.5,0.375,55),ncol=1)

  lpZ<-simZ(1e6)%*%matrix(c(0.69,0.69,0.69,0.35,-0.035),ncol=1)

  x<-0:30
  lpX0<-x*delta0
  for (i in 1:length(chpnts)){
    if (i==1) { lpX0=lpX0+delta1* ifelse(x-chpnts[i]<0,0, x-chpnts[i])
    } else {
      if (alter==TRUE&(i/2==floor(i/2))) lpX0=lpX0-delta1* ifelse(x-chpnts[i]<0,0, x-chpnts[i])   else lpX0=lpX0+delta1* ifelse(x-chpnts[i]<0,0, x-chpnts[i])
    }
  }

  lpX0<-rep(lpX0,1e6/31)
  lp<-as.vector(lpZ)+lpX0


  beta0.opt<-uniroot(get.beta0_pois, c(-200,200), tol=1e-12, lp=lp, target.mean=mn,maxiter=1e6)$root #for p=1 it would be enough to do it once, since we get the same throughout, but for p>1 we get different optimum constant
  beta0.opt
}

#logistic and beta

get.beta0_logit<-function(beta0, lp, target.prop){
  -target.prop + mean(1/ (1+exp(-beta0 - lp)))
}

get_intercept_logistic<-function(mn,delta0,delta1,chpnts,alter){
  #lpZ<-c(0.69,0.69,0.69,0.35,-0.035)%*%matrix(c(0.158,0.158,0.5,0.375,55),ncol=1)

  lpZ<-simZ(1e6)%*%matrix(c(0.69,0.69,0.69,0.35,-0.035),ncol=1)

  x<-0:30
  lpX0<-x*delta0
  for (i in 1:length(chpnts)){
    if (i==1) { lpX0=lpX0+delta1* ifelse(x-chpnts[i]<0,0, x-chpnts[i])
    } else {
      if (alter==TRUE&(i/2==floor(i/2))) lpX0=lpX0-delta1* ifelse(x-chpnts[i]<0,0, x-chpnts[i])  else lpX0=lpX0+delta1* ifelse(x-chpnts[i]<0,0, x-chpnts[i])
    }
  }

  lpX0<-rep(lpX0,1e6/31)
  lp<-as.vector(lpZ)+lpX0

  #lp<-as.vector(lpZ+lpX0)

  beta0.opt<-uniroot(get.beta0_logit, c(-200,200), tol=1e-8, lp=lp, target.prop=mn)$root #for p=1 it would be enough to do it once, since we get the same throughout, but for p>1 we get different optimum constant
  beta0.opt
}


###functions to get output

est.coef<-function(fit){

  list<-rep(NA,12)
  names(list) <- c("(Intercept)", " ZZ1", " ZZ2"," ZZ3"," ZZ4"," ZZ5","x","U1.x","U2.x","U3.x","U4.x","U5.x")
  results<-list
  h<- ifelse(is.null(nrow(fit$psi)),0,nrow(fit$psi))
  coef<-coef(fit)
  i<-1
  while(i<h+8){
    results[i]<-coef[i]
    i=i+1}
  return(results)
}

est.nochp<-function(fit){
  h<- ifelse(is.null(nrow(fit$psi)),0,nrow(fit$psi))
  return(h)
}

est.location<-function(fit){
  t <-NA
  h<- ifelse(is.null(nrow(fit$psi)),0,nrow(fit$psi))
  if (h==0) return (t)
  else return (fit$psi[,"Est."])
}

est.location2<-function(x){
  xx<-rep(NA,5)
  xx[1:length(x)]<-x
  xx
}

est.CI<-function(fit,tci){ #tci can only be "gradient" since we are also considering nonlinear models
  ci<-matrix(NA,ncol=2,nrow=5)
  h<- ifelse(is.null(nrow(fit$psi)),0,nrow(fit$psi))
  if (h==0) return (ci)
  if (class(fit)[1]!="rq"){
    if(tci=="delta") ci1<-confint.segmented(fit,method="delta")
    if(tci=="score") ci1<-confint.segmented(fit,method="score")
    if (tci=="gradient") ci1<-confint.segmented(fit,method="gradient")
  } else {
    if(tci=="delta") ci1<-confint.segmented(fit,method="delta",.vcov =vcov.rq(fit))
    if(tci=="score") ci1<-confint.segmented(fit,method="score",.vcov =vcov.rq(fit))
    if (tci=="gradient") ci1<-confint.segmented(fit,method="gradient",.vcov =vcov.rq(fit))
  }
  ci[1:nrow(ci1),]<-ci1[,2:3]
  return(ci)

}

my.chp.true<-function(x){
  xx<-rep(NA,5)
  xx[1:length(x)]<-x
  xx
}







####################
#delta0=0.2
#delta1=0
#task="medium"
#no=ii
#alter=kk
#n=500
#model="linear"
#var.sigma=3
#phi.beta=NULL
#intercept=1
#mean_y=.50
#max_segments=5

###only makes sense for linear! weird results for either AIC/BIC for QR and beta!

##only implemented for linear model!
#i=m
#delta0=0.2
#delta1=0.6
#task="medium"
#alter=kk
#no=ii
#n=250
#model=model

#var.sigma=3
#phi.beta=phi.beta
#intercept=intercept
#mean_y=mean_y

F_simulated<-function(M,data,npsi1){
  y<-data$Y

  fit<-lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=data)

  n<-nrow(data)

  npsi2<-npsi1+1


  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  fit.seg2<-try(segmented(fit, seg.Z=~X, npsi=npsi2,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)


  if (class(fit.seg)[1]=="try-error") {
    pv0<-NA

    } else {


   if (class(fit.seg2)[1]=="try-error") {

     pv0<-1
   }  else {


      X<-data$X
      Z<-model.matrix(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=data)

      if (npsi1!=0&!is.null(fit.seg$psi)){
        if(nrow(fit.seg$psi)!=0){
          for (jj in 1:nrow(fit.seg$psi)){

            Z<-cbind( Z, ifelse(X-fit.seg$psi[jj,2]<0,0, X-fit.seg$psi[jj,2]),-ifelse(X>fit.seg$psi[jj,2],1,0))

          }
        }
      }

      Z2<-model.matrix(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=data)

      if (npsi2!=0&!is.null(fit.seg2$psi)){
        if(nrow(fit.seg2$psi)!=0){
          for (jj in 1:nrow(fit.seg2$psi)){

            Z2<-cbind( Z2, ifelse(X-fit.seg2$psi[jj,2]<0,0, X-fit.seg2$psi[jj,2]),-ifelse(X>fit.seg2$psi[jj,2],1,0))

          }
        }
      }

      pp<-try(solve(t(Z)%*%Z),silent=TRUE)
      pp2<-try(solve(t(Z2)%*%Z2),silent=TRUE)

      #fit.seg


      if (class(pp)[1]=="try-error"|class(pp2)[1]=="try-error") {

        pv0=1



      } else {
        I<-matrix(0,n,n)
        diag(I)<-1

        P0<-Z%*%pp%*%t(Z)
        P1<-Z2%*%pp2%*%t(Z2)

        res<-(I-P0)%*%matrix(y,ncol=1)



        F<- y%*%(I-P0)%*%matrix(y,ncol=1)/y%*%(I-P1)%*%matrix(y,ncol=1)

        F<-c(F)
        if (F==1) pv0=1 else {

          Fi<-list()



          for (ii  in 1:M){

            resids<-matrix( res*rnorm(n)  ,ncol=1)

            num<-  t(resids)%*%(I-P0)%*%resids
            denom<-(y%*%P0+ t(resids)  )%*%(I-P1)%*%(P0%*%matrix(y,ncol=1)+resids)

            Fi[[ii]]<-c(num/denom)





          } #end for ii

          pvF<-(sum( as.numeric(Fi)>=F)+1)/(length(Fi)+1)

          pv0<-pvF
          #pv0<-mm*pv0#if we adjust the p
        } #end if F

      } #end if clas fit.seg
    } #end if class
  }
  return(pv0)

}




gof_simulated<-function(M,data,npsi1){
  df<-data


  fit<-lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=data)


  if (npsi1==0) fit.seg<-fit else   fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,
                                                           control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)


  if (class(fit.seg)[1]=="try-error") {
    pv0<-NA


  } else {


      #sbst<-which(x>fit.seg$psi[nrow(fit.seg$psi),2]) # not very promissing, it is ok now and much better then the overall test
      res<-fit.seg$residuals
      res.sd<-res/sd(res)
      res.sd.sim<-res.sd

      dfo<-data.frame(res.sd,X=data$X)
      setDT(dfo)
      dfo[, res.sd := mean(res.sd, na.rm = TRUE), by = X]

      res.sd<-dfo$res.sd


      sbst<-order(data$X)



      W<-1/sqrt(length(res.sd))*cumsum(res.sd[sbst])



      #plot((1:n)[sbst],W,type="l",ylim=c(-1,1))
      #abline(v=fit.seg$psi[,2])
      Wi<-list()
      X<-df$X
      Z<-model.matrix(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X, data=df)
      if (!is.null(fit.seg$psi)){
        for (jj in 1:nrow(fit.seg$psi)){

          Z<-cbind( Z, ifelse(X-fit.seg$psi[jj,2]<0,0, X-fit.seg$psi[jj,2]),-ifelse(X>fit.seg$psi[jj,2],1,0))
          #Z<-cbind( Z,  (x-fit.seg$psi[jj,2])*ifelse( x>fit.seg$psi[jj,2],1,0  )  )

        }
      }
      #Z2<-cbind(Z,-ifelse(x>fit.seg$psi[1,2],1,0))
      #pp<-solve(t(Z)%*%Z)
      pp<-try(solve(t(Z)%*%Z),silent=TRUE) #this makes all the difference! now it is not general, can only be a single changepoint! I think the design matrix for more changepoints is then different!
      #  pp<-try( (t(Z)%*%Z)^(-1),silent=TRUE) #this makes all the difference! now it is not general, can only be a single changepoint! I think the design matrix for more changepoints is then different!

      if (class(pp)[1]=="try-error") {
        pv0<-1


      } else {

        for (ii in 1:M){


          #     res<-fit.seg$res



          rr<-res.sd.sim*rnorm(length(res.sd))

          #rr<-res.sd*(2*rbinom(n,size=1,prob=0.5)  -1)
          dfor<-data.frame(rr,X=data$X)
          setDT(dfor)
          dfor[, rr := mean(rr, na.rm = TRUE), by = X]

          rr<-dfor$rr


          #rm<-matrix(rr,ncol=ncol(Z),nrow=length(rr))

          prvi.del<-cumsum(rr[sbst])


          drugi.del.1<-apply(Z[sbst, ],2,cumsum)

          drugi.del.2<-pp%*%t(Z[sbst, ])%*%matrix( rr[sbst],ncol=1)



          #drugi.del.1<-apply(Z2[sbst, ],2,cumsum)

          #drugi.del.2<-pp%*%t(Z2[sbst, ])%*%matrix( rr[sbst],ncol=1)




          drugi.del<-drugi.del.1%*%drugi.del.2



          Wi[[ii]]<-1/sqrt(length(prvi.del))*( prvi.del-drugi.del)
          #lines((1:n)[sbst],Wi[[ii]],col="grey")




        } #end for ii


        # lines((1:n)[sbst],W,col="red")
        CvM<-function(x){ sum( x**2  )  }
        ts<-CvM(W)
        tsi<-unlist(lapply(Wi,CvM))

        pv<-(sum(tsi>=ts,na.rm=T)+1)/(sum(!is.na(tsi),na.rm=T)+1)
        #mtext(paste("p=",round(pv,6)))

        pv0<-pv
        #pv0<-(mm+1)*pv0#if we adjust the p
      } #end ifelse 2
     #end ifelse 1
  }
  return(pv0)

}



#######################main function to perform simulation For now just on one fix data set

my_sim<-function(i,delta0,
                 delta1,
                 task,
                 no,
                 alter,
                 n,
                 model,
                 var.sigma,phi.beta,
                 intercept,mean_y,max_segments=5,data=mydata){



    int=intercept


  #mydata<-mysim_data(n=n,dif=task,no=no, delta0=delta0,delta1=delta1,alter=alter,model=model,intercept=int,var.sigma=var.sigma,phi.beta = phi.beta)
  mydata<-data
  # hist(mydata$data$Y)
  #boxplot(mydata$data$Y~mydata$data$X)
  mnY<-mean(mydata$data$Y)
  sdY<-sd(mydata$data$Y)

    fit<-try(lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=mydata$data),silent=TRUE)


  if (class(fit)[1]!="try-error"){
    conv.init=1




    ##GoF and F based procedure

    pv.gof.sim <-  pv.f.sim <-  rep(NA, max_segments)

    time.gof.sim  <-  time.f.sim <-  rep(NA, max_segments)

    for (npsi in 0:(max_segments - 1)) {

      t <- system.time({
        pv.gof.sim[npsi + 1] <- gof_simulated(
          M = 1000, data = mydata$data, npsi1 = npsi
        )
      })
      time.gof.sim[npsi + 1] <- t["elapsed"]



      t <- system.time({
        pv.f.sim[npsi + 1] <- F_simulated(
          M = 1000, data = mydata$data, npsi1 = npsi
        )
      })
      time.f.sim[npsi + 1] <- t["elapsed"]



    }


    get_num<-function(pv,alpha){
      if (is.na(pv[1])) NA else {
      if (pv[1]>alpha|is.na(pv[2])) 0 else {
        if (pv[2]>alpha|is.na(pv[3])) 1 else {
          if (pv[3]>alpha|is.na(pv[4])) 2 else {
            if (pv[4]>alpha|is.na(pv[5])) 3 else {
              if (pv[5]>alpha) 4 else 5
            }
          }
        }
      }
      }
    }

    n_gof_sim<-get_num(pv.gof.sim,0.05)
     n_f_sim<-get_num(pv.f.sim,0.05)

    n_gof_sim_bonf<-get_num(pv.gof.sim,0.05/max_segments)
     n_f_sim_bonf<-get_num(pv.f.sim,0.05/max_segments)




  } else {
    conv.init=0


    n_gof_sim<-NA
     n_f_sim<-NA

    n_gof_sim_bonf<-NA
     n_f_sim_bonf<-NA

    pv.gof.sim <-   pv.f.sim <-   rep(NA, max_segments)

    time.gof.sim  <-   time.f.sim <-   rep(NA, max_segments)



  }



  desc.data<-c(mnY,sdY,no,my.chp.true(mydata$chpoints),conv.init,int)
  names(desc.data)<-c("mean_Y","sd_Y","true_n_chpnt",paste("true_chpnt",1:5,sep="_"),"convergence_init","true_intercept")

  result<-c(
            n_gof_sim,

            n_f_sim,


            n_gof_sim_bonf,

            n_f_sim_bonf
            )
  names(result)<-c("gofsim","fsim",
                   "gofsimbonf","fsimbonf")
  #result


  nm<-paste("d0",delta0,"d1",delta1,task,"a",ifelse(alter==TRUE,1,0),"no",no,
            "n",n,model,sep=".")

    nm<-paste(nm,"vs",var.sigma,"in",intercept,sep=".")


  nm<-paste(nm,".txt",sep="")

  nmdes<-paste("results/d",nm,sep=".")
  nmr1<-paste("results/r1",nm,sep=".")

  write(desc.data,append=T,ncol=length(desc.data),sep="\t",file=nmdes)

  write(result,append=T,ncol=length(result),sep="\t",file=nmr1)

  #save p-values and times!
res_p<-c(pv.gof.sim ,   pv.f.sim  )
time_p<-c(time.gof.sim  ,   time.f.sim  )

nmr2<-paste("results/pval",nm,sep=".")
nmr3<-paste("results/times",nm,sep=".")

write(res_p,append=T,ncol=length(res_p),sep="\t",file=nmr2)
write(time_p,append=T,ncol=length(time_p),sep="\t",file=nmr3)


}






