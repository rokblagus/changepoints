
###here I tried to check why F-P performs so poorly



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
    eps<-rnorm(n,sd=sqrt(var.sigma))
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








F_permuted_setpsi<-function(M,data,npsi1,chpoints){

  fit<-lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=data)
  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,psi=chpoints,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  fit.seg1<-try(segmented(fit, seg.Z=~X, npsi=npsi1+1,psi=c(chpoints,15),   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)

  if (class(fit.seg)[1]=="try-error") return(rep(NA,4)) else {
    if(class(fit.seg1)[1]=="try-error") return(c(rep(1,3),NA)) else {

      fited<-predict(fit.seg)
      res<-data$Y-fited
      #res<-res-mean(res)
      #fited<-fit.seg$fitted.values
      # res<-fit.seg$residuals
      res1<-data$Y-predict(fit.seg1)
      #res1<-res1-mean(res1)
      #res1<-fit.seg1$residuals

      n<-length(res)
      F<-sum(res**2)/sum(res1**2)

      Fi<-rep(NA,M)

      datai<-data
      seeds<-1:20000


      for (i in 1:M){
        if (i>20000)  set.seed(max(seeds)+i-20000)  else      set.seed(seeds[i])

        yr<-fited+res[sample(1:n)]

        datai$Y<-yr

        fiti<- lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=datai)
        if (npsi1==0) fit.segi<-fiti else fit.segi<-try(segmented(fiti, seg.Z=~X, npsi=npsi1,psi=chpoints,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
        #dont refit H0 model
        #fit.segi<-fit.seg

        fit.seg1i<-try(segmented(fiti, seg.Z=~X, npsi=npsi1+1,psi=c(chpoints,15),   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)

        if (class(fit.segi)[1]=="try-error"|class(fit.seg1i)[1]=="try-error"){

          Fi[i]<-F

        } else {
          resi<-datai$Y-predict(fit.segi)
          res1i<-datai$Y-predict(fit.seg1i)
          #resi<-resi-mean(resi)
          #res1i<-res1i-mean(res1i)
          #resi<-fit.segi$residuals
          #res1i<-fit.seg1i$residuals

          Fi[i]<-sum(resi**2)/sum(res1i**2)

        }

      }
      #hist(Fi)
      #abline(v=F)
      getp<-function(x,y){  y<-c(x,y); mean(y>=x,na.rm=TRUE) }
      getp2<-function(x,y){  y<-c(x,y); mean(y>x,na.rm=TRUE) }
      getp3<-function(x,y){  y<-c(x,y[y!=x]); mean(y>=x,na.rm=TRUE) }
      return(c(getp(F,Fi),getp2(F,Fi),getp3(F,Fi),sum(Fi==F,na.rm=TRUE)))

    }
  }
}


gof_permuted_setpsi<-function(M,data,npsi1,chpoints){

  fit<-lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=data)
  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,psi=chpoints,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  if (class(fit.seg)[1]=="try-error") return(rep(NA,4)) else {
    fited<-predict(fit.seg)
    res.perm<-data$Y-fited
    #res<-res-mean(res)
    #fited<-fit.seg$fitted.values
    res<-fit.seg$residuals
    res.sd<-res/sd(res)

    #library(data.table)
    df<-data.frame(res.sd,X=data$X)
    setDT(df)
    df[, res.sd := mean(res.sd, na.rm = TRUE), by = X]


    or<-order(data$X)
    n<-length(res)
    W<-cumsum(df$res.sd[or])
    Wi<-list()

    datai<-data
    #i=0
    #while (i < M){
    # i=i+1
    seeds<-1:20000

    for (i in 1:M){
      if (i>20000)  set.seed(max(seeds)+i-20000)  else      set.seed(seeds[i])
      yr<-fited+res.perm[sample(1:n)]

      datai$Y<-yr

      fiti<- lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=datai)

      if (npsi1==0) fit.segi<-fiti else fit.segi<-try(segmented(fiti, seg.Z=~X, npsi=npsi1,psi=chpoints,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
      if (class(fit.segi)[1]!="try-error"){
        #resi<-datai$Y-predict(fit.segi)
        #resi<-resi-mean(resi)
        resi<-fit.segi$residuals


        res.sdi<-resi/sd(resi)

        dfi<-data.frame(res.sdi,X=data$X)
        setDT(dfi)
        dfi[, res.sdi := mean(res.sdi, na.rm = TRUE), by = X]



        Wi[[i]]<-cumsum(dfi$res.sdi[or])
      } else {Wi[[i]]<-W}

    }
    cvm<-function(x) sum(x**2)
    ts<-cvm(W)
    tsi<-unlist(lapply(Wi,cvm))
    #hist(tsi)
    #abline(v=ts)
    getp<-function(x,y){  y<-c(x,y); mean(y>=x,na.rm=TRUE) }
    getp2<-function(x,y){  y<-c(x,y); mean(y>x,na.rm=TRUE) }
    getp3<-function(x,y){  y<-c(x,y[y!=x]); mean(y>=x,na.rm=TRUE) }
    return(c(getp(ts,tsi),getp2(ts,tsi),getp3(ts,tsi),sum(tsi==ts,na.rm=TRUE)))

  }
}

F_permuted<-function(M,data,npsi1){

  fit<-lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=data)
  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  fit.seg1<-try(segmented(fit, seg.Z=~X, npsi=npsi1+1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)

  if (class(fit.seg)[1]=="try-error") return(rep(NA,4)) else {
    if(class(fit.seg1)[1]=="try-error") return(c(rep(1,3),NA)) else {

      fited<-predict(fit.seg)
      res<-data$Y-fited
      #res<-res-mean(res)
      #fited<-fit.seg$fitted.values
      # res<-fit.seg$residuals
      res1<-data$Y-predict(fit.seg1)
      #res1<-res1-mean(res1)
      #res1<-fit.seg1$residuals

      n<-length(res)
      F<-sum(res**2)/sum(res1**2)

      Fi<-rep(NA,M)

      datai<-data
      seeds<-1:20000


      for (i in 1:M){
        if (i>20000)  set.seed(max(seeds)+i-20000)  else      set.seed(seeds[i])

        yr<-fited+res[sample(1:n)]

        datai$Y<-yr

        fiti<- lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=datai)
        if (npsi1==0) fit.segi<-fiti else fit.segi<-try(segmented(fiti, seg.Z=~X, npsi=npsi1,  model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
        #dont refit H0 model
        #fit.segi<-fit.seg

        fit.seg1i<-try(segmented(fiti, seg.Z=~X, npsi=npsi1+1,  model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)

        if (class(fit.segi)[1]=="try-error"|class(fit.seg1i)[1]=="try-error"){

          Fi[i]<-F

        } else {
          resi<-datai$Y-predict(fit.segi)
          res1i<-datai$Y-predict(fit.seg1i)
          #resi<-resi-mean(resi)
          #res1i<-res1i-mean(res1i)
          #resi<-fit.segi$residuals
          #res1i<-fit.seg1i$residuals

          Fi[i]<-sum(resi**2)/sum(res1i**2)

        }

      }
      #hist(Fi)
      #abline(v=F)
      getp<-function(x,y){  y<-c(x,y); mean(y>=x,na.rm=TRUE) }
      getp2<-function(x,y){  y<-c(x,y); mean(y>x,na.rm=TRUE) }
      getp3<-function(x,y){  y<-c(x,y[y!=x]); mean(y>=x,na.rm=TRUE) }
      return(c(getp(F,Fi),getp2(F,Fi),getp3(F,Fi),sum(Fi==F,na.rm=TRUE)))



    }
  }
}


gof_permuted<-function(M,data,npsi1){

  fit<-lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=data)
  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  if (class(fit.seg)[1]=="try-error") return(rep(NA,4)) else {
    fited<-predict(fit.seg)
    res.perm<-data$Y-fited
    #res<-res-mean(res)
    #fited<-fit.seg$fitted.values
    res<-fit.seg$residuals
    res.sd<-res/sd(res)

    #library(data.table)
    df<-data.frame(res.sd,X=data$X)
    setDT(df)
    df[, res.sd := mean(res.sd, na.rm = TRUE), by = X]


    or<-order(data$X)
    n<-length(res)
    W<-cumsum(df$res.sd[or])
    Wi<-list()

    datai<-data
    #i=0
    #while (i < M){
    # i=i+1
    seeds<-1:20000

    for (i in 1:M){
      if (i>20000)  set.seed(max(seeds)+i-20000)  else      set.seed(seeds[i])
      yr<-fited+res.perm[sample(1:n)]

      datai$Y<-yr

      fiti<- lm(Y~Z.Z1+Z.Z2+Z.Z3+Z.Z4+Z.Z5+X,data=datai)

      if (npsi1==0) fit.segi<-fiti else fit.segi<-try(segmented(fiti, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
      if (class(fit.segi)[1]!="try-error"){
        #resi<-datai$Y-predict(fit.segi)
        #resi<-resi-mean(resi)
        resi<-fit.segi$residuals


        res.sdi<-resi/sd(resi)

        dfi<-data.frame(res.sdi,X=data$X)
        setDT(dfi)
        dfi[, res.sdi := mean(res.sdi, na.rm = TRUE), by = X]



        Wi[[i]]<-cumsum(dfi$res.sdi[or])
      } else {Wi[[i]]<-W}

    }
    cvm<-function(x) sum(x**2)
    ts<-cvm(W)
    tsi<-unlist(lapply(Wi,cvm))
    #hist(tsi)
    #abline(v=ts)
    getp<-function(x,y){  y<-c(x,y); mean(y>=x,na.rm=TRUE) }
    getp2<-function(x,y){  y<-c(x,y); mean(y>x,na.rm=TRUE) }
    getp3<-function(x,y){  y<-c(x,y[y!=x]); mean(y>=x,na.rm=TRUE) }
    return(c(getp(ts,tsi),getp2(ts,tsi),getp3(ts,tsi),sum(tsi==ts,na.rm=TRUE)))

  }
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

        npsi=no


      pv.f.perm.set <- try(F_permuted_setpsi(
        M = 1000, data = mydata$data, npsi1 = npsi,chpoints=mydata$chpoints
      ),silent=TRUE)

      pv.f.perm <- try(F_permuted(
        M = 1000, data = mydata$data, npsi1 = npsi
      ),silent=TRUE)

      pv.gof.perm<-try(gof_permuted(M=1000,data=mydata$data,npsi1= npsi),silent=TRUE)

      pv.gof.perm.set<-try(gof_permuted_setpsi(M=1000,data=mydata$data,npsi1= npsi,chpoints=mydata$chpoints),silent=TRUE)

      if (class(pv.f.perm.set)[1]=="try-error") pv.f.perm.set<-rep(NA,4)

      if (class(pv.f.perm)[1]=="try-error") pv.f.perm<-rep(NA,4)
      if (class(pv.gof.perm)[1]=="try-error") pv.gof.perm<-rep(NA,4)
      if (class(pv.gof.perm.set)[1]=="try-error") pv.gof.perm.set<-rep(NA,4)

      res.cfs<-c(pv.f.perm.set,pv.f.perm,pv.gof.perm,pv.gof.perm.set)


      nm<-paste("d0",delta0,"d1",delta1,task,"a",ifelse(alter==TRUE,1,0),"no",no,
                "n",n,model,sep=".")

      nm<-paste(nm,"vs",var.sigma,"in",intercept,sep=".")


      nm<-paste(nm,".txt",sep="")

      nmr1<-paste("results/r1",nm,sep=".")

      write(res.cfs,append=T,ncol=length(res.cfs),sep="\t",file=nmr1)





}



