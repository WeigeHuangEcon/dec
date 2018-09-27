#' @title obd
#' @description do Oaxaca-Blinder decomposition
#' @param formula a formula: y ~ x1+x2
#' @param ref reference group. ref=g1 or g2.
#' @param data1 data for group 1
#' @param data2 data for group 2
#' @return aggregate and detailed decomposition components
#' @examples
#' data(five)
#' na=five[,1:6]
#' eu=five[,7:12]
#' colnames=c("Mkt.RF","SMB","HML","RMW", "CMA","SMALL.LoBM")
#' colnames(na)=colnames
#' colnames(eu)=colnames
#' formula=SMALL.LoBM~Mkt.RF+SMB+HML+RMW+CMA
#' res=obd(formula,ref="g1",data1=na,data2=eu)
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

#' @title sd_obd
#' @description compute standard deviation of Oaxaca-Blinder decomposition
#' @inheritParams obd
#' @param B the number of bootstrap iterations
#' @return standard deviations of aggregate and detailed decomposition components
#' @examples
#' data(five)
#' na=five[,1:6]
#' eu=five[,7:12]
#' colnames=c("Mkt.RF","SMB","HML","RMW", "CMA","SMALL.LoBM")
#' colnames(na)=colnames
#' colnames(eu)=colnames
#' formula=SMALL.LoBM~Mkt.RF+SMB+HML+RMW+CMA
#' B=10
#' sd=sd_obd(formula,ref="g1",data1=na,data2=eu,B=B)
#' @export
#'
#'
sd_obd=function(formula,ref="g1",data1,data2,B){
  drr=pbapply::pbsapply(1:B, function(i) {
    data1=data1[sample(nrow(data1),replace = T),]
    data2=data2[sample(nrow(data2),replace = T),]
    obd(formula,ref=ref,data1,data2)
  })
  sd=apply(drr, 1, sd)
  sd
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
#' @examples
#' data(five)
#' na=five[,1:6]
#' colnames=c("Mkt.RF","SMB","HML","RMW", "CMA","SMALL.LoBM")
#' colnames(na)=colnames
#' formula=SMALL.LoBM~Mkt.RF+SMB+HML+RMW+CMA
#' y=na$SMALL.LoBM
#' y0=seq(quantile(y,probs = 0.05),quantile(y,probs = 0.95),length=100)
#' model=dr(formula,y0=y0,data=na,link="logit")
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
#' @examples
#' data(five)
#' na=five[,1:6]
#' eu=five[,7:12]
#' colnames=c("Mkt.RF","SMB","HML","RMW", "CMA","SMALL.LoBM")
#' colnames(na)=colnames
#' colnames(eu)=colnames
#' formula=SMALL.LoBM~Mkt.RF+SMB+HML+RMW+CMA
#' y=na$SMALL.LoBM
#' y0=seq(quantile(y,probs = 0.05),quantile(y,probs = 0.95),length=100)
#' cfs=cfs_dr(formula,y0,ref="g1",data1=na,data2=eu,link="logit")
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
#' @examples
#' data(five)
#' na=five[,1:6]
#' eu=five[,7:12]
#' colnames=c("Mkt.RF","SMB","HML","RMW", "CMA","SMALL.LoBM")
#' colnames(na)=colnames
#' colnames(eu)=colnames
#' formula=SMALL.LoBM~Mkt.RF+SMB+HML+RMW+CMA
#' y=na$SMALL.LoBM
#' y0=seq(quantile(y,probs = 0.05),quantile(y,probs = 0.95),length=100)
#' cfs=cfs_dr(formula,y0,ref="g1",data1=na,data2=eu,link="logit")
#' qs=seq(0.1,0.9,length.out = 10)
#' qsres=qs_dr(object=cfs,y0=y0,qs=qs)
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
#' @examples
#' data(five)
#' na=five[,1:6]
#' eu=five[,7:12]
#' colnames=c("Mkt.RF","SMB","HML","RMW", "CMA","SMALL.LoBM")
#' colnames(na)=colnames
#' colnames(eu)=colnames
#' formula=SMALL.LoBM~Mkt.RF+SMB+HML+RMW+CMA
#' y=na$SMALL.LoBM
#' y0=seq(quantile(y,probs = 0.05),quantile(y,probs = 0.95),length=100)
#' qs=seq(0.1,0.9,length.out = 10)
#' deres=de_dr(formula,y0,ref="g1",data1=na,data2=eu,link="logit",qs=qs)
#'
#' @export
#'
de_dr=function(formula,y0,ref="g1",data1,data2,link="logit",qs) {
  cfsr_dr=cfs_dr(formula,y0,ref=ref,data1,data2,link=link)
  qsr_dr=qs_dr(object=cfsr_dr,y0,qs)
  dr=dec_dr(object=qsr_dr)
  formula=as.formula(formula)
  dat1=model.frame(terms(formula,data=data1),data=data1)
  Overall=rowSums(dr)
  Agg_c=rowSums(dr[,1:(ncol(dat1)-1)])
  Agg_s=Overall-Agg_c
  drs=cbind(Overall,Agg_c,Agg_s,dr)
  colnames(drs)=c("Overall","Agg_c","Agg_s",colnames(dr))
  drs
}

#' @title sd_dr
#' @description compute standard deviation for components
#' @inheritParams de_dr
#' @param B the number of bootstrap iterations
#' @return standard deviations for decomposition components
#' @examples
#' data(five)
#' na=five[,1:6]
#' eu=five[,7:12]
#' colnames=c("Mkt.RF","SMB","HML","RMW", "CMA","SMALL.LoBM")
#' colnames(na)=colnames
#' colnames(eu)=colnames
#' formula=SMALL.LoBM~Mkt.RF+SMB+HML+RMW+CMA
#' y=na$SMALL.LoBM
#' y0=seq(quantile(y,probs = 0.05),quantile(y,probs = 0.95),length=100)
#' qs=seq(0.1,0.9,length.out = 10)
#' B=3
#' sd_dr(formula,y0,ref="g1",data1=na,data2=eu,link="logit",qs,B)
#' @export
#'
#'

sd_dr=function(formula,y0,ref="g1",data1,data2,link="logit",qs,B){
  drr=pbapply::pblapply(1:B, function(i) {
    data1=data1[sample(nrow(data1),replace = T),]
    data2=data2[sample(nrow(data2),replace = T),]
    de_dr(formula,y0,ref="g1",data1,data2,link="logit",qs)
  })
  sd=apply(simplify2array(drr), 1:2, sd)
  sd
}

#' @title de_dr.plot
#' @description plot the reselts from de_dr
#' @param object decomposition components from de_dr
#' @param qs quantiles to compute
#' @param kind kinds of decomposition.It could be "agg", "composition", "structure".
#' @export
de_dr.plot=function(object,kind="agg",qs){
  if(kind=="agg"){
    matplot(qs,object[,1:3],type = "l",lwd = 3,xlab = "Quantile",ylab = "Effects")
    legend("topleft",c("Overall","Composition","Structure"),col = 1:3,lty = 1:3,lwd=3,cex = 0.5)
  }
  if(kind=="composition"){
    n=(ncol(object)-2)/2 +2
    cc=colnames(object)[4:n]
    dd=object[,c(2,4:n)]
    matplot(qs,dd,type = "l",lwd = 3,xlab = "Quantile",ylab = "Effects")
    legend("topleft",c("Agg_c",cc),col = 1:length(dd),lty = 1:length(dd),lwd=3,cex = 0.5)
  }
  if(kind=="structure"){
    n=(ncol(object)-2)/2 +2
    cc=colnames(object)[-(1:n)]
    ccc=c("Agg_s",cc)
    dd=cbind(object[,3],object[,-(1:n)])
    matplot(qs,dd,type = "l",lwd = 3,xlab = "Quantile",ylab = "Effects")
    legend("topleft",ccc,col = 1:length(dd),lty = 1:length(dd),lwd=3,cex = 0.5)
  }
}


