# =====================================================================
# Authors: Ana Radulovic and Rok Blagus
# Part of the code for: "Estimating the number of changepoints in segmented
# regression: new sequential tests and comparative evaluation"
#
# functions_tests.R - F-type and GoF tests (simulation- and permutation-based,
# with and without the adjustment for heteroscedasticity), adapted from the
# functions used in the simulation study (final_gof_F_functions.R).
#
# Changes with respect to the simulation code (the logic of the tests and the
# handling of failed fits are unchanged):
#  1. model formula Y ~ SEX + AGE_YRS + X instead of Y ~ Z.Z1 + ... + Z.Z5 + X;
#  2. F_simulated, F_simulated_type2: computed without n x n projection
#     matrices (P0 v = Z (Z'Z)^{-1} Z'v), i.e. the same algebra;
#  3. heteroscedasticity-adjusted versions: the Gaussian kernel on X is
#     computed over the unique values of X instead of an n x n matrix
#     (identical result);
#  4. bandwidth: the simulation code uses 450/sqrt(n), where 450 = 30^2/2 for
#     X on [0, 30]. Here BW_CONST/sqrt(n) with BW_CONST = (range of X)^2/2
#     (1989-2021: 32^2/2 = 512), multiplied by BW_FACTOR (1 by default;
#     0.5 and 2 in the sensitivity analysis). BW_CONST <- 450 reproduces the
#     simulation code;
#  5. the "example of use" at the end of the original file was removed.
#
# Within the permutation-based tests, a different seed is set for every
# permutation; a seed should also be set before every call of a test.
# =====================================================================

####required packages

library(segmented)
library(data.table)

#this line forces data.table to use a single core

setDTthreads(1)



####settings for the real data

BW_CONST<-NULL   # NULL = (range of X)^2/2 from the data; 450 = simulation code
BW_FACTOR<-1     # 0.5 or 2 for the sensitivity analysis

bw_real<-function(data){
  bc<-if (is.null(BW_CONST)) diff(range(data$X))^2/2 else BW_CONST
  BW_FACTOR*bc/sqrt(nrow(data))
}

#local residual variance: sigma.u_i = sum_j K(x_i,x_j) res_j^2 / sum_j K(x_i,x_j)
#K(a,b) = exp(-(a-b)^2/bandwidth), as gausskernel(); computed over the unique values of x
sigma_u_by_x<-function(x,res,bandwidth){
  ux<-sort(unique(x))
  g<-match(x,ux)
  S<-tapply(res^2,factor(g,levels=seq_along(ux)),sum)
  cnt<-tabulate(g,nbins=length(ux))
  Kx<-exp(-outer(ux,ux,"-")^2/bandwidth)
  s<-as.vector(Kx%*%S)/as.vector(Kx%*%cnt)
  s[g]
}

#####aux functions

#needed for hetero-adjusted methods
gausskernel<-function (X = NULL, sigma = NULL){
  return(exp(-1 * as.matrix(dist(X)^2)/sigma))
}

est.nochp<-function(fit){
  h<- ifelse(is.null(nrow(fit$psi)),0,nrow(fit$psi))
  return(h)
}




####F-P


F_permuted<-function(M,data,npsi1){

  fit<-lm(Y~SEX+AGE_YRS+X,data=data)
  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  fit.seg1<-try(segmented(fit, seg.Z=~X, npsi=npsi1+1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)

  if (class(fit.seg)[1]=="try-error") NA else {
    if(class(fit.seg1)[1]=="try-error") 1 else {

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

        fiti<- lm(Y~SEX+AGE_YRS+X,data=datai)
        if (npsi1==0) fit.segi<-fiti else fit.segi<-try(segmented(fiti, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
        #dont refit H0 model
        #fit.segi<-fit.seg

        fit.seg1i<-try(segmented(fiti, seg.Z=~X, npsi=npsi1+1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)

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
      return(getp(F,Fi))

    }
  }
}

####F-P heteroscedasticity adjusted


F_permuted_hetero<-function(M,data,npsi1){

  bandwidth<-bw_real(data)

  fit<-lm(Y~SEX+AGE_YRS+X,data=data)
  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  fit.seg1<-try(segmented(fit, seg.Z=~X, npsi=npsi1+1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)

  if (class(fit.seg)[1]=="try-error") NA else {
    if(class(fit.seg1)[1]=="try-error") 1 else {

      fited<-predict(fit.seg)
      res<-data$Y-fited
      #res<-res-mean(res)
      #fited<-fit.seg$fitted.values
      # res<-fit.seg$residuals
      res1<-data$Y-predict(fit.seg1)
      #res1<-res1-mean(res1)
      #res1<-fit.seg1$residuals

      n<-length(res)

      sigma.u<-sigma_u_by_x(data$X,res,bandwidth)  # = n x n kernel of the simulation code, computed by year


      F<-sum(res**2)/sum(res1**2)

      Fi<-rep(NA,M)

      datai<-data
      seeds<-1:20000


      for (i in 1:M){
        if (i>20000)  set.seed(max(seeds)+i-20000)  else      set.seed(seeds[i])

        yr<-fited+((res/sqrt(sigma.u))[sample(1:n)])*sqrt(sigma.u)

        datai$Y<-yr

        fiti<- lm(Y~SEX+AGE_YRS+X,data=datai)
        if (npsi1==0) fit.segi<-fiti else fit.segi<-try(segmented(fiti, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
        #dont refit H0 model
        #fit.segi<-fit.seg

        fit.seg1i<-try(segmented(fiti, seg.Z=~X, npsi=npsi1+1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)

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
      return(getp(F,Fi))

    }
  }
}



####F-S

F_simulated<-function(M,data,npsi1){
  y<-data$Y

  fit<-lm(Y~SEX+AGE_YRS+X,data=data)

  n<-nrow(data)

  npsi2<-npsi1+1


  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  fit.seg2<-try(segmented(fit, seg.Z=~X, npsi=npsi2,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)


  if (class(fit.seg)[1]=="try-error"|class(fit.seg2)[1]=="try-error") {
    pv0<-NA


  } else {
    if (est.nochp(fit.seg2)!=npsi2 |est.nochp(fit.seg)!=npsi1) {
      pv0<-1
    } else {
      X<-data$X
      Z<-model.matrix(Y~SEX+AGE_YRS+X,data=data)

      if (npsi1!=0&!is.null(fit.seg$psi)){
        if(nrow(fit.seg$psi)!=0){
          for (jj in 1:nrow(fit.seg$psi)){

            Z<-cbind( Z, ifelse(X-fit.seg$psi[jj,2]<0,0, X-fit.seg$psi[jj,2]),-ifelse(X>fit.seg$psi[jj,2],1,0))

          }
        }
      }

      Z2<-model.matrix(Y~SEX+AGE_YRS+X,data=data)

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

        pv0=NA



      } else {
        # without n x n matrices: P0 %*% v = Z %*% (pp %*% (t(Z) %*% v)); same algebra as in the simulation code
        P0v<-function(v) Z%*%(pp%*%crossprod(Z,v))
        P1v<-function(v) Z2%*%(pp2%*%crossprod(Z2,v))
        ym<-matrix(y,ncol=1)
        P0y<-P0v(ym)

        res<-ym-P0y



        F<- sum(ym*(ym-P0y))/sum(ym*(ym-P1v(ym)))

        F<-c(F)
        if (F==1) pv0=1 else {

          Fi<-list()



          for (ii  in 1:M){

            resids<-matrix( res*rnorm(n)  ,ncol=1)

            num<-  sum(resids*(resids-P0v(resids)))
            ww<-P0y+resids
            denom<-sum(ww*(ww-P1v(ww)))

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


####F-S2


F_simulated_type2<-function(M,data,npsi1){
  y<-data$Y

  fit<-lm(Y~SEX+AGE_YRS+X,data=data)

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
      Z<-model.matrix(Y~SEX+AGE_YRS+X,data=data)

      if (npsi1!=0&!is.null(fit.seg$psi)){
        if(nrow(fit.seg$psi)!=0){
          for (jj in 1:nrow(fit.seg$psi)){

            Z<-cbind( Z, ifelse(X-fit.seg$psi[jj,2]<0,0, X-fit.seg$psi[jj,2]),-ifelse(X>fit.seg$psi[jj,2],1,0))

          }
        }
      }

      Z2<-model.matrix(Y~SEX+AGE_YRS+X,data=data)

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
        # without n x n matrices: P0 %*% v = Z %*% (pp %*% (t(Z) %*% v)); same algebra as in the simulation code
        P0v<-function(v) Z%*%(pp%*%crossprod(Z,v))
        P1v<-function(v) Z2%*%(pp2%*%crossprod(Z2,v))
        ym<-matrix(y,ncol=1)
        P0y<-P0v(ym)

        res<-ym-P0y



        F<- sum(ym*(ym-P0y))/sum(ym*(ym-P1v(ym)))

        F<-c(F)
        if (F==1) pv0=1 else {

          Fi<-list()



          for (ii  in 1:M){

            resids<-matrix( res*rnorm(n)  ,ncol=1)

            num<-  sum(resids*(resids-P0v(resids)))
            ww<-P0y+resids
            denom<-sum(ww*(ww-P1v(ww)))

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




###############


################



###Gof-S


gof_simulated<-function(M,data,npsi1){
  df<-data


  fit<-lm(Y~SEX+AGE_YRS+X,data=data)


  if (npsi1==0) fit.seg<-fit else   fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,
                                                           control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)


  if (class(fit.seg)[1]=="try-error") {
    pv0<-NA


  } else {

    if (est.nochp(fit.seg)!=npsi1){
      pv0<-1
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
      Z<-model.matrix(Y~SEX+AGE_YRS+X, data=df)
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
    } #end ifelse 1
  }
  return(pv0)

}


###Gof-P


gof_permuted<-function(M,data,npsi1){

  fit<-lm(Y~SEX+AGE_YRS+X,data=data)
  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  if (class(fit.seg)[1]=="try-error") NA else {
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

      fiti<- lm(Y~SEX+AGE_YRS+X,data=datai)

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
    return(getp(ts,tsi))

  }
}


###Gof-P heteroscedasticity adjusted


gof_permuted_hetero<-function(M,data,npsi1){
  bandwidth<-bw_real(data)

  fit<-lm(Y~SEX+AGE_YRS+X,data=data)
  if (npsi1==0) fit.seg<-fit else fit.seg<-try(segmented(fit, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
  if (class(fit.seg)[1]=="try-error") NA else {
    fited<-predict(fit.seg)
    res.perm<-data$Y-fited
    #res<-res-mean(res)
    #fited<-fit.seg$fitted.values
    res<-fit.seg$residuals


    sigma.u<-sigma_u_by_x(data$X,res,bandwidth)  # = n x n kernel of the simulation code, computed by year


    res.sd<-res/sqrt(sigma.u)

    #library(data.table)
    df<-data.frame(res.sd,X=data$X)
    setDT(df)
    df[, res.sd := mean(res.sd, na.rm = TRUE), by = X]


    or<-order(data$X)
    n<-length(res)
    W<-cumsum(df$res.sd[or])/sqrt(n)
    Wi<-list()

    datai<-data
    #i=0
    #while (i < M){
    # i=i+1
    seeds<-1:20000

    for (i in 1:M){
      if (i>20000)  set.seed(max(seeds)+i-20000)  else      set.seed(seeds[i])
      yr<-fited+((res.perm/sqrt(sigma.u))[sample(1:n)])*sqrt(sigma.u)

      datai$Y<-yr

      fiti<- lm(Y~SEX+AGE_YRS+X,data=datai)

      if (npsi1==0) fit.segi<-fiti else fit.segi<-try(segmented(fiti, seg.Z=~X, npsi=npsi1,   model = TRUE,control=seg.control(n.boot = 50, alpha = 0.01) ),silent=TRUE)
      if (class(fit.segi)[1]!="try-error"){
        #resi<-datai$Y-predict(fit.segi)
        #resi<-resi-mean(resi)
        resi<-fit.segi$residuals
        fitedi<-predict(fit.segi)
        sigma.ui<-sigma_u_by_x(data$X,resi,bandwidth)  # = n x n kernel of the simulation code, computed by year


        res.sdi<-resi/sqrt(sigma.ui)

        dfi<-data.frame(res.sdi,X=data$X)
        setDT(dfi)
        dfi[, res.sdi := mean(res.sdi, na.rm = TRUE), by = X]



        Wi[[i]]<-cumsum(dfi$res.sdi[or])/sqrt(n)
      } else {Wi[[i]]<-W}

    }
    cvm<-function(x) sum(x**2)
    ts<-cvm(W)
    tsi<-unlist(lapply(Wi,cvm))
    #hist(tsi)
    #abline(v=ts)
    getp<-function(x,y){  y<-c(x,y); mean(y>=x,na.rm=TRUE) }
    return(getp(ts,tsi))

  }
}

###GoF-S2



gof_simulated_type2<-function(M,data,npsi1){
  df<-data


  fit<-lm(Y~SEX+AGE_YRS+X,data=data)


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
    Z<-model.matrix(Y~SEX+AGE_YRS+X, data=df)
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


######################
