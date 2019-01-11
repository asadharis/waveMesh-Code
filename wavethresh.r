wavethresh = function(y, thresh.value="universal", sigma=NA, p=1, low.level=4, by.level=F, filter.number=2, family="DaubExPhase", bc="periodic"){
    # Wavelet-based smoothing by Lp penalized least squares
    # INPUT
    # y = noisy signal
    # thresh.value = "universal", "minimax", "SURE", "BIC", "SLIC"
    # sigma = standard deviation of noise
    # p in [0,1] (0=hard, 1=soft)
    # low.level = lowest level where to start thresholding
    # by.level = T (levelwise) or F (global)


    if(!is.na(p)){
      if(p==1) type="soft"
      if(p==0) type="hard"
    }
    else type="lpin01"

    N=length(y)
    if(!is.numeric(sigma)){
      ytHaar=wd(y,1,family="DaubExPhase")  
      sigma=mad(accessD(ytHaar,log(N,2)-1),center=0)
    }
    y=y/sigma
    yt=wd(y, filter.number=filter.number, family=family, bc=bc)

    if(by.level==F){
      if(thresh.value=="universal") thresh.value=sqrt(2*log(N))
      if(thresh.value=="minimax") {
        Ns=2^c(6:16)
        tmmS=c(1.474,1.669,1.859,2.045,2.226,2.403,2.575,2.743,2.906,3.066,3.221)
        tmmH=c(2.697,2.913,3.117,3.312,3.497,3.674,3.844,4.008,4.166,4.319,4.467)
        if(type=="soft") thresh.value=tmmS[match(N,Ns)]
        if(type=="hard") thresh.value=tmmH[match(N,Ns)]
      }
      if((thresh.value=="SURE")&(type=="soft")){
        aMLE=rev(yt$D)
        ## DJ 1995 JASA's Hybrid method eq. 14 page 1208
        if(mean(aMLE^2-1)>(log(N,2)^(3/2)/sqrt(N))){
          aMLE=aMLE[-c(1:(2^low.level-1))]
          aMLEabssort = sort(abs(aMLE))
          aMLEabssort2 = aMLEabssort^2
          lambdas=aMLEabssort;
          ND=length(aMLE)
          tp1 = cumsum(aMLEabssort2)+((ND-1):0)*aMLEabssort2;
          risk_hatS = (ND - (2 * (1:ND)) + tp1)/ND;
          thresh.value = lambdas[which.min(risk_hatS)]
        }
        else thresh.value=sqrt(2*log(N))
      }
      if((thresh.value=="SLIC")&(type=="soft")){
        ND=length(yt$D)
        lambdaU=sqrt(2*log(ND))
        lambda0=lambdaU/2
        a.hat=lp_shrink(yt$D, 1, lambda0)
        l1mean=mean(abs(a.hat))
        again=T
        while(again){
          lambda1=lambdaU/(1+l1mean*lambdaU)
          again=((abs(lambda1-lambda0)/lambda1)>1.e-6)
          lambda0=lambda1
          a.hat=lp_shrink(yt$D, 1, lambda0)
          l1mean=mean(abs(a.hat))
        }
        thresh.value=lambda0
      }
      if((thresh.value=="BIC")&(type=="hard")) thresh.value=sqrt(log(N))

      ytthresh=threshold.wd(yt, levels=low.level:(yt$nlevels-1), type=type, policy="manual", value=thresh.value)
      ytthresh$D=ytthresh$D*sigma
      ytthresh$C=ytthresh$C*sigma
    }
    else{
      ###########################  
      ## Levelwise thresholding
      ytthresh=yt
      for(level in low.level:(yt$nlevel-1)){
        ytj=accessD(yt,level)
        Nj=length(ytj)
        if(thresh.value=="universal") thresh.value.j=sqrt(2*log(Nj))
        if(thresh.value=="minimax") {
          Ns=2^c(3:16)
          # The first 4 values are extrapolated ones.
          tmmS=c(0.674,0.874,1.074,1.274,1.474,1.669,1.859,2.045,2.226,2.403,2.575,2.743,2.906,3.066,3.221)
          tmmH=c(1.757,2.007,2.247,2.477,2.697,2.913,3.117,3.312,3.497,3.674,3.844,4.008,4.166,4.319,4.467)
          if(type=="soft") thresh.value.j=tmmS[match(Nj,Ns)]
          if(type=="hard") thresh.value.j=tmmH[match(Nj,Ns)]
        }
        if((thresh.value=="SURE")&(type=="soft")){
          ## DJ 1995 JASA's Hybrid method eq. 14 page 1208
          if(mean(ytj^2-1)>(log(Nj,2)^(3/2)/sqrt(Nj))){
            aMLE=ytj
            aMLEabssort = sort(abs(aMLE))
            aMLEabssort2 = aMLEabssort^2
            lambdas=aMLEabssort;
            tp1 = cumsum(aMLEabssort2)+((Nj-1):0)*aMLEabssort2;
            risk_hatS = (Nj - (2 * (1:Nj)) + tp1)/Nj;
            thresh.value.j = lambdas[which.min(risk_hatS)]
          }
          else thresh.value.j=sqrt(2*log(Nj))
        }
        if((thresh.value=="SURE")&(type=="hard")){
          ## DJ 1995 JASA's Hybrid method eq. 14 page 1208
          if(mean(ytj^2-1)>(log(Nj,2)^(3/2)/sqrt(Nj))){
            aMLE=ytj
            aMLEabssort = sort(abs(aMLE))
            aMLEabssort2 = aMLEabssort^2
            risk_hatH=(Nj-(2 * (1:Nj)) + cumsum(aMLEabssort2))/Nj;
            ahatU=ytj*(abs(ytj)>sqrt(2*log(Nj)));
            for(n in 1:Nj){
              phi=aMLEabssort[n];
              risk_hatH[n]=risk_hatH[n]+2*phi/sqrt(2*pi)*mean(exp(-0.5*(-phi-ahatU)^2)+exp(-0.5*(phi-ahatU)^2))
            }
            thresh.value.j = aMLEabssort[which.min(risk_hatH)]
          }
          else thresh.value.j=sqrt(2*log(Nj))
        }
        if((thresh.value=="SLIC")&(type=="soft")){
          lambdaU=sqrt(2*log(Nj))
          lambda0=lambdaU/2
          ytj.hat=lp_shrink(ytj, 1, lambda0)
          l1mean=mean(abs(ytj.hat))
          again=T
          while(again){
            lambda1=lambdaU/(1+l1mean*lambdaU)
            again=((abs(lambda1-lambda0)/lambda1)>1.e-6)
            lambda0=lambda1
            ytj.hat=lp_shrink(ytj, 1, lambda0)
            l1mean=mean(abs(ytj.hat))
          }
          thresh.value.j=lambda0
        }
        if((thresh.value=="BIC")&(type=="hard")) thresh.value.j=sqrt(log(Nj))
        if((thresh.value=="SURE")&is.na(p)){
           outmin=nlminb(.5, SURE_lp, control=list(iter.max=20, x.tol=.05), lower=0, upper=1, yt=ytj, out_nb=1)
           p.j=outmin$par
           out=SURE_lp(p.j, ytj, 2)
           lambda.j=out$lambda_SURE
           thresh.value.j=(2-p.j)*(lambda.j*(2*(1-p.j))^(p.j-1))^(1/(2-p.j))
        }
        if((thresh.value=="SLIC")&is.na(p)){
          p1=1/log(Nj);
          phiU=sqrt(2*log(Nj));
          lambdapU=(phiU/(2-p1))^(2-p1)/(2-2*p1)^(p1-1)
          dyl = (2-2*p1)^((p1-1)/(2-p1))*lambdapU^((p1-1)/(2-p1));
          dyll= (2-2*p1)^((p1-1)/(2-p1))*(p1-1)/(2-p1)*lambdapU^((2*p1-3)/(2-p1));
          taupU=p1*phiU*lambdapU*(dyl)^2/(Nj*dyl+p1*lambdapU*dyll);
          ct0=phiU/taupU;ct1=0;low=lambdapU;high=lambdapU
          while(rootlambda(low,p1,ct0,ct1,Nj)<=0) {low=low/2}
          while(rootlambda(high,p1,ct0,ct1,Nj)>=0) {high=high*2}  
          lambda1=uniroot(rootlambda,p=p1,ct0=ct0,ct1=ct1,N=Nj,interval=c(low,high));lambda1=lambda1$root  
  
          MAXITER=20;
          again_p=1; again_lambda=1; iter=0;
          while((again_p|again_lambda)&(iter<MAXITER)) {
            iter=iter+1;
            lambda0=lambda1;
            p0=p1;
            ytj.hat = lp_shrink(ytj, p0, lambda0);
            outmin=nlminb(.5,condloglikinp,lower=1.e-6,upper=1,alpha=ytj.hat,lambda=lambda0,N=Nj,tau=taupU);
            p1=outmin$par
            again_p=(mean(abs((p1-p0)/p1))>.001);
        
            lambdapU=(phiU/(2-p1))^(2-p1)/(2-2*p1)^(p1-1); 
            dyl=(2-2*p1)^((p1-1)/(2-p1))*lambdapU^((p1-1)/(2-p1));
            dyll=(2-2*p1)^((p1-1)/(2-p1))*(p1-1)/(2-p1)*lambdapU^((2*p1-3)/(2-p1));        
            taupU=p1*phiU*lambdapU*(dyl)^2/(Nj*dyl+p1*lambdapU*dyll);
        
            Lp=mean((abs(ytj.hat))^p1);
            ct0=phiU/taupU;ct1=Lp;
            low=lambda0;
            while(rootlambda(low,p1,ct0,ct1,Nj)<=0) {low=low/2}
            high=lambdapU;
            while(rootlambda(high,p1,ct0,ct1,Nj)>=0) {high=high*2}
            lambda1=uniroot(rootlambda,p=p1,ct0=ct0,ct1=ct1,N=Nj,interval=c(low,high));lambda1=lambda1$root  
        
            again_lambda=(mean(abs((lambda1-lambda0)/lambda1))>.01);
          }        
          lambda.j=lambda1;
          p.j=p1;
          thresh.value.j=(2-p.j)*(lambda.j*(2*(1-p.j))^(p.j-1))^(1/(2-p.j))
        }
        if(!is.na(p)&((p==0)||(p==1))) ytthresh=threshold.wd(ytthresh, levels=level, type=type, policy="manual", value=thresh.value.j)
        else {
          ytthresh.j=lp_shrink(ytj, p.j,lambda.j)
          ytthresh=putD(ytthresh, level, ytthresh.j)
        }
      }
      ytthresh$D=ytthresh$D*sigma
      ytthresh$C=ytthresh$C*sigma
    }

    mu.hat=wr(ytthresh)

    out=NULL
    out$alpha.MLE=yt
    out$alpha.hat=ytthresh
    out$mu.hat=mu.hat
    out$sigma.hat=sigma
    return(out)
}

##########################################
## Two functions employed by SLIC: rootlambda and condloglikinp
rootlambda=function(lambda,p,ct0,ct1,N){

    tau=sqrt(2*log(N))/ct0;
    if(lambda<0){
      fl=1.e10;
    }
    else {
        dyll=(2-2*p)^((p-1)/(2-p))*(p-1)/(2-p)*lambda^((2*p-3)/(2-p));
        temp1=(2-2*p)^((p-1)/(2-p))*lambda^((p-1)/(2-p));
        fl=-ct1*lambda + 1/p - ct0*lambda/N*temp1 + lambda/N*dyll/temp1;
    }
    return(fl)
    }
    
###############################################    
condloglikinp=function(p, alpha,lambda, N, tau){

  J=log(length(alpha))/log(2);
  Nm=length(alpha);
  Lp=sum((abs(alpha))^p);
  fp=Nm*lgamma(1/p+1)-Nm/p*log(lambda)+lambda*Lp;
  phi=(2-p)*(lambda*(2-2*p)^(p-1))^(1/(2-p));
  
  phiN=sqrt(2*log(N));
  dN=phiN-(log(log(N))+log(4*pi)-2*log(2))/2/phiN;
  ct1=phiN*(phi/tau-dN);
  logG1=log(phiN/tau)-exp(-ct1)-ct1;
  dyl=(2-2*p)^((p-1)/(2-p))*lambda^((p-1)/(2-p))
  logG1=logG1+log(dyl);

  fp=fp-logG1;    

  return(fp)
}
