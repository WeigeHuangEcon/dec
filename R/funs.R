#' @title obd
#' @description do Oaxaca-Blinder decomposition
#' @param formula a formula: y ~ x1+x2
#' @param ref reference group. ref=g1 or g2.
#' @param data1 data for group 1
#' @param data2 data for group 2
#' @return aggregate and detailed decomposition components
#' @export
#'
#'

obd=function(formula,ref="g1",data1,data2){
  formula=as.formula(formula)
  dat1=model.frame(terms(formula,data=data1),data=data1)
  dat2=model.frame(terms(formula,data=data2),data=data2)
  model1=glm(formula,data = data1)
  model2=glm(formula,data = data2)
  coef1=model1$coefficients
  coef2=model2$coefficients
  factors=names(model1$coefficients)
  datt1=dat1[,c(factors[-1])]
  datt2=dat2[,c(factors[-1])]

  if(ref=="g1") {
    dc=rep(0,length(factors))
    for (i in 2:length(factors)) {
      dc[i]=(colMeans(datt1)-colMeans(datt2))[i-1]*coef1[i]
    }
    dc=dc[-1]
    ds=rep(0,length(factors))
    ds[1]=coef1[1]-coef2[1]
    for (i in 2:length(factors)) {
      ds[i]=colMeans(datt2)[i-1]*(coef1[i]-coef2[i])
    }
   c=sum(dc)
   s=sum(ds)
   d=c+s
   dr=c(d,c,s,dc,ds)
   cnames=paste0(factors[-1],sep="_c")
   snames=paste0(c("Constant",factors[-1]),sep = "_s")
   names(dr)=c("Overall","Agg_c","Agg_s",cnames,snames)
   dr
  }
  if(ref=="g2") {
    dc=rep(0,length(factors))
    for (i in 2:length(factors)) {
      dc[i]=(colMeans(datt1)-colMeans(datt2))[i-1]*coef2[i]
    }
    dc=dc[-1]
    ds=rep(0,length(factors))
    ds[1]=coef1[1]-coef2[1]
    for (i in 2:length(factors)) {
      ds[i]=colMeans(datt1)[i-1]*(coef1[i]-coef2[i])
    }
    c=sum(dc)
    s=sum(ds)
    d=c+s
    dr=c(d,c,s,dc,ds)
    cnames=paste0(factors[-1],sep="_c")
    snames=paste0(c("Constant",factors[-1]),sep = "_s")
    names(dr)=c("Overall","Agg_c","Agg_s",cnames,snames)
    dr
  }
  dr
}

#' @title ind
#' @description create indicator matrix
#' @param y outcome variable
#' @param y0 a grid of threshold values
#' @return indicator matrix
#' @export
#'
#'
ind=function(y,y0){
  indicatormatrix=matrix(0,nrow = length(y),ncol = length(y0))
  for(i in 1:length(y0)){
    indicatormatrix[,i]=1*(y<=y0[i])
  }
  indicatormatrix
}

#' @title dr
#' @description distribution regression
#' @param formula a formula: y ~ x1+x2
#' @param y0 a grid of threshold values
#' @param data data
#' @param link link function: logit or probit
#' @return a list of models for y0
#' @export
#'
#'
dr=function(formula,y0,data,link="logit"){
  formula=as.formula(formula)
  dat=model.frame(terms(formula,data=data),data=data)
  y=dat[,1]
  indicatormatrix=ind(y,y0)
  model=list(list())
  for(i in 1:length(y0)){
    ff=paste(colnames(dat)[-1],collapse = "+")
    formulaa=as.formula(paste("indicatormatrix[,i]~",ff))
    model[[i]]=glm(formulaa,family = binomial(link=link),data=dat)
  }
  model
}


#' @title cfs_dr
#' @description compute (counterfactual) distributions using distribution regression
#' @param formula a formula: y ~ x1+x2
#' @param y0 a grid of threshold values
#' @param ref reference group. ref=g1 or g2.
#' @param data1 data for group 1
#' @param data2 data for group 2
#' @param link link function: logit or probit
#' @return a series of (counterfactual) distributions
#' @export
#'
#'
cfs_dr=function(formula,y0,ref="g1",data1,data2,link="logit"){
  formula=as.formula(formula)
  dat1=model.frame(terms(formula,data=data1),data=data1)
  dat2=model.frame(terms(formula,data=data2),data=data2)
  x1=cbind(1,dat1[,-1])
  x2=cbind(1,dat2[,-1])
  rs_dr1=dr(formula,y0,data1,link=link)
  rs_dr2=dr(formula,y0,data2,link=link)
  if(ref=="g1"){
    if(link=="logit"){
      logitF=function(x){
        logitF=exp(x)/(1+exp(x))
        logitF
      }
      fs=matrix(0,nrow = length(y0),ncol =2*ncol(x1) )
      for (i in 1:length(y0)) {
        coef1=rs_dr1[[i]]$coefficients
        coef2=rs_dr2[[i]]$coefficients
        #f1=mean(logitF(as.matrix(x1) %*% coef1))
        #f2=mean(logitF(as.matrix(x2) %*% coef2))
        f_x12=rep(0,ncol(x1))
        for (j in 1:ncol(x1)) {
          x12j=cbind(x2[,1:j],x1[,-c(1:j)])
          f_x12[j]=mean(logitF(as.matrix(x12j) %*% coef1),na.rm = T)
        }
        f_coe12=rep(0,ncol(x1))
        for (j in 1:ncol(x1)) {
          coe12j=c(coef2[1:j],coef1[-(1:j)])
          f_coe12[j]=mean(logitF(as.matrix(x2) %*% coe12j),na.rm = T)
        }
        fs[i,]=c(f_x12,f_coe12)
      }
      fs
    }
    if(link=="probit"){
      fs=matrix(0,nrow = length(y0),ncol =2*ncol(x1) )
      for (i in 1:length(y0)) {
        coef1=rs_dr1[[i]]$coefficients
        coef2=rs_dr2[[i]]$coefficients
        #f1=mean(logitF(as.matrix(x1) %*% coef1))
        #f2=mean(logitF(as.matrix(x2) %*% coef2))
        f_x12=rep(0,ncol(x1))
        for (j in 1:ncol(x1)) {
          x12j=cbind(x2[,1:j],x1[,-c(1:j)])
          f_x12[j]=mean(pnorm(as.matrix(x12j) %*% coef1),na.rm = T)
        }
        f_coe12=rep(0,ncol(x1))
        for (j in 1:ncol(x1)) {
          coe12j=c(coef2[1:j],coef1[-(1:j)])
          f_coe12[j]=mean(pnorm(as.matrix(x2) %*% coe12j),na.rm = T)
        }
        fs[i,]=c(f_x12,f_coe12)
      }
     fs
    }
  }

  if(ref=="g2"){
    if(link=="logit"){
      logitF=function(x){
        logitF=exp(x)/(1+exp(x))
        logitF
      }
      fs=matrix(0,nrow = length(y0),ncol =2*ncol(x1) )
      for (i in 1:length(y0)) {
        coef1=rs_dr1[[i]]$coefficients
        coef2=rs_dr2[[i]]$coefficients
        #f1=mean(logitF(as.matrix(x1) %*% coef1))
        #f2=mean(logitF(as.matrix(x2) %*% coef2))
        f_x12=rep(0,ncol(x1))
        for (j in 1:ncol(x1)) {
          x12j=cbind(x2[,1:j],x1[,-c(1:j)])
          f_x12[j]=mean(logitF(as.matrix(x12j) %*% coef1),na.rm = T)
        }
        f_coe12=rep(0,ncol(x1))
        for (j in 1:ncol(x1)) {
          coe12j=c(coef2[1:j],coef1[-(1:j)])
          f_coe12[j]=mean(logitF(as.matrix(x2) %*% coe12j),na.rm = T)
        }
        fs[i,]=c(f_x12,f_coe12)
      }
      fs
    }
    if(link=="probit"){
      fs=matrix(0,nrow = length(y0),ncol =2*ncol(x1) )
      for (i in 1:length(y0)) {
        coef1=rs_dr1[[i]]$coefficients
        coef2=rs_dr2[[i]]$coefficients
        #f1=mean(logitF(as.matrix(x1) %*% coef1))
        #f2=mean(logitF(as.matrix(x2) %*% coef2))
        f_x12=rep(0,ncol(x1))
        for (j in 1:ncol(x1)) {
          x12j=cbind(x2[,1:j],x1[,-c(1:j)])
          f_x12[j]=mean(pnorm(as.matrix(x12j) %*% coef2),na.rm = T)
        }
        f_coe12=rep(0,ncol(x1))
        for (j in 1:ncol(x1)) {
          coe12j=c(coef2[1:j],coef1[-(1:j)])
          f_coe12[j]=mean(pnorm(as.matrix(x1) %*% coe12j),na.rm = T)
        }
        fs[i,]=c(f_x12,f_coe12)
      }
      fs
    }
  }
  cnames=paste0(colnames(dat1)[-1],sep="_c")
  snames=paste0(c("Constant",colnames(dat1)[-1]),sep="_s")
  colnames(fs)=c("f0",cnames,snames)
 fs
}



#' @title qs_dr
#' @description inverse the distribution function to obtain quantiles
#' @param object results from cfs_dr
#' @param y0 a grid of threshold values
#' @param qs quantiles to compute
#' @return a series of quantiles corresponding to each (counterfactual) distribution
#' @export
#'

qs_dr=function(object,y0,qs){
  f2q=function(fs,y0,qs){
    qss=rep(0,length(qs))
    for (i in 1:length(qs)) {
      qss[i]= y0[which(sort(fs) >= qs[i])[1]]
    }
    qss
  }
  qs_dr=pbapply::pbsapply(1:ncol(object),function(i) f2q(object[,i],y0,qs=qs))
  colnames(qs_dr)=colnames(object)
  qs_dr
}


#' @title dec_dr
#' @description decompose quantile differences
#' @param object results from qs_dr
#' @return decomposition components of the quantile differences
#' @keywords internal
#' @export
#'

dec_dr=function(object){
  d=matrix(0,nrow = nrow(object),ncol = ncol(object)-1)
  for (i in 1:ncol(d)) {
    d[,i]=object[,i]-object[,i+1]
  }
  colnames(d)=colnames(object)[-1]
  d
}

#' @title de_dr
#' @description decompose quantile differences using distribution regression
#' @inheritParams cfs_dr
#' @inheritParams dr
#' @inheritParams qs_dr
#' @return decomposition components of the quantile differences
#' @export
#'
de_dr=function(formula,y0,ref="g1",data1,data2,link="logit",qs) {
  cfsr_dr=cfs_dr(formula,y0,ref=ref,data1,data2,link=link)
  qsr_dr=qs_dr(object=cfsr_dr,y0,qs)
  dr=dec_dr(object=qsr_dr)
  dr
}


