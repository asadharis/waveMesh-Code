AMlet = function(y,X,thresh.value="universal", sigma=NA, p=1, low.level=4, by.level=F, filter.number=2, family="DaubExPhase", bc="periodic", conv.thresh=1.e-10, max.iter=ncol(X)*1000){
  ## Fits additive models using wavethresh

    N=length(y)
    P=ncol(X)
    continue=T
    Xo=apply(X,2,order)
    Xr=apply(X,2,rank)

    ## linear model fit
    outlm=lm(y~X, model=T)
    alpha0=outlm$coef[1]
    Shat=t(outlm$coef[-1]*t(X))
    res=y-alpha0-apply(Shat,1,sum)
    iter=0
    sigmahat=mad(outlm$res)
    if(is.numeric(sigma)) sigmahat=sigma
    while(continue&(iter<max.iter)){
      iter=iter+1
      Shat0=Shat
      sigmahatnew=c()
      for(pp in 1:P){
        resp=res+Shat[,pp]
        respo=resp[Xo[,pp]]
        outp=wavethresh(respo, thresh.value, sigmahat, p, low.level, by.level, filter.number, family, bc)
        Shatp=outp$mu.hat[Xr[,pp]]
        res=resp-Shatp
        Shat[,pp]=Shatp
        if(!is.numeric(sigma)){
          ytHaar=wd(respo,1,family="DaubExPhase")  
          sigmahatnew=c(sigmahatnew,mad(accessD(ytHaar,log(N,2)-1),center=0))
        }
      }
      if(!is.numeric(sigma)) sigmahat=(sigmahat+median(sigmahatnew))/2
      else sigmahat=sigma
      alpha0new=mean(res+alpha0)
      res=res+alpha0-alpha0new
      alpha0=alpha0new
      continue=((mean(abs(Shat-Shat0))/mean(abs(Shat)))>conv.thresh)
    }

  out=NULL
  out$Shat=Shat
  out$alpha0=alpha0
  out$muhat=alpha0+apply(Shat,1,sum)
  out$sigmahat=sigmahat
  out$Xo=Xo
  out$outlm=outlm
  out$converged=(iter<max.iter)
  return(out)
}
