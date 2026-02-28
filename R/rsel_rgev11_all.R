library(goftest)
library(Rsolnp)
library(lmomco)
library(eva)

#--------------------------------------------------------------------
rsel.rgev11= function(xdat, model='rgev11', sigL=0.05, 
                          method=c("ed","ccdf","spacing"),
                          num_inits=15, seq.cut=TRUE, 
                          qqplot=FALSE){
  
  z=list()                                     # model='rgev11' 
#  work.mle=paste(model,".fit.park",sep="")     # work.mle='rgev11.fit.park'
  z$model= model
  z$method= method
  
  ns=nrow(xdat); dim=ncol(xdat)
  marg.pval= rep(NA,dim)
  mle=matrix(NA,dim,5)
  iD=mtheoU=matrix(NA, ns,dim)
  result = matrix(NA,dim,8)
  
  if(seq.cut==TRUE){                   # seq.cut= FALSE
    
    sth=0
    while (sth < dim){                 # sth=1, sth=2  # qqplot=TRUE
      sth= sth+1 
      
      if(sth==1) {st.para=NULL
      }else{st.para= result[sth-1,2:6]}
      result[sth,1]= sth
      
#      cat("sth=",sth,"\n")
      
      gof = rgev11.fit.park(xdat, r= sth, 
                            num_inits=num_inits,
                            start.para=st.para)
      result[sth,2:6]= gof$mle
      
      if(result[sth,6] < -0.5){
        gof = rgev11.fit.park(xdat, r= sth, 
                              num_inits=num_inits*2,
                              start.para=st.para, penk="CD")
        result[sth,2:6]= gof$mle
      }
      
      result[sth,8]= gof$nllh
      
      ztilda = ns.trsf.rlos.Gum(xdat,para=result[sth,2:6],
                                sth=sth, method=method)
     
      if(method=="ed" & sth ==1) result[sth,7]=1
      if(method=="ed" & sth >=2){
        
        Ed.test=ed.gum.gof(ztilda,r=sth,parx=c(0,1,-1e-5))
        result[sth,7] = Ed.test$p.value
        
        if(result[sth,7] < sigL){z$r.sel=r.sel= sth-1; break}
        
      }

      if(method=="ccdf"){
        
        mtheoU[,sth]= rcond.gof.gev11(xdat, r=sth, model=model,
                                      parx=result[sth,2:6])
        
        # H0: r= sth,  H1: r= sth-1
        result[sth,7]= cvm.test(mtheoU[,sth], null="punif")$p.value
        
         #       cat("sth, p.val=", sth, result[sth,6],"\n" )
        
        if(result[sth,7] < sigL){z$r.sel=r.sel= sth-1; break}
      }
      
      if(method=="spacing" & sth ==1) result[sth,7]=1
      if(method=="spacing" & sth >= 2){
        
        iD[,sth]= spacing.gof.gev11(xdat, r= sth-1, 
                                    parx=result[sth,2:6])
        
        result[sth,7]= cvm.test(iD[,sth], null="pexp")$p.value

        #       cat("sth, p.val=", sth, result[sth,6],"\n" )

        if(result[sth,7] < sigL){z$r.sel=r.sel= sth-1; break}
      }

    } #end while
    
  }else if(seq.cut==FALSE){
    
    for(sth in 1:dim){
      
      if(sth==1) {st.para=NULL
      }else{st.para= result[sth-1,2:6]}
      result[sth,1]= sth
      
#      cat("sth=",sth,"\n")
      
      gof = rgev11.fit.park(xdat, r= sth, 
                            num_inits=num_inits,
                            start.para=st.para)
      result[sth,2:6]= gof$mle
      
      if(result[sth,6] < -0.5){
        gof = rgev11.fit.park(xdat, r= sth, 
                              num_inits=num_inits*2,
                              start.para=st.para, penk="CD")
        result[sth,2:6]= gof$mle
      }
      
      result[sth,8]= gof$nllh
      
      ztilda = ns.trsf.rlos.Gum(xdat,para=result[sth,2:6],
                                sth=sth, method=method)
      
      if(method=="ed" & sth >=2){
        
        Ed.test=ed.gum.gof(ztilda,r=sth,parx=c(0,1,-1e-5))
        result[sth,7] = Ed.test$p.value

      }

      if(method=="ccdf"){
        
        mtheoU[,sth]= rcond.gof.gev11(xdat, r=sth, model=model,
                                      parx=result[sth,2:6])
        result[sth,7]= cvm.test(mtheoU[,sth], null="punif")$p.value
        
      }
      if(method=="spacing" & sth >= 2){
        
        iD[,sth]= spacing.gof.gev11(xdat, r= sth-1, 
                                    parx=result[sth,2:6])
        result[sth,7]= cvm.test(iD[,sth], null="pexp")$p.value 
      }
      
    } #end for sth
    
    if(method=="spacing" | method=="ed") result[1,7]=1
    id.rej = which(result[,7] < sigL)
    z$r.sel = min(id.rej, dim+1)-1
    
  } #end if seq.cut
  
  if(seq.cut==TRUE){ work.r = min(z$r.sel+1,dim)
  }else{ work.r = dim}
  
  if(qqplot==TRUE){
    ff=(seq(1:ns)-0.35)/ns
    par(new=F); par(mfrow=c(2,3))

    if(method=="ccdf"){
      for (sth in 1:work.r){
        plot(ff, mtheoU[,sth], main=paste("r=",sth,sep=" "),
             ylab="ccdf.model", 
             xlab=paste("EDF: p=",round(result[sth,7],3),sep=" "),
             pch=16) 
        abline(a=0,b=1) }
    }
    if(method=="spacing"){
      for(sth in 2:work.r){
        plot(ff, pexp(iD[,sth]), main=paste("r=",sth,sep=" "),
             ylab="spacing.model", 
             xlab=paste("EDF: p=",round(result[sth,7],3),sep=" "),
             pch=16, col="blue") 
        abline(a=0,b=1) }
    }
  } # end if qqplot
  
  colnames(result) = c("r","est.loc0","est.loc1","est.sc0",
                       "est.sc1","est.shape",
                       "p.value","nllh")
  
  if(is.null(z$r.sel)) z$r.sel= dim
  z$result = na.omit(as.data.frame(result))
  z$source= 'rsel.gev11'
  z
}

#-----------------------
ns.trsf.rlos.Gum=function(xdat, para=null, sth=null,
                          method=NULL)
{
  ns=nrow(xdat)
  if(is.null(maxr)) maxr=ncol(xdat)
  mu0 =para[1]
  mu1 =para[2]
  sig0=para[3]
  sig1=para[4]
  xi  =para[5]   # hosking style
  
  mut <- mu0+mu1*c(1:ns)
  sigt<- exp(sig0+sig1*c(1:ns))
  xi  <- xi 
  ztilda=matrix(0,ns,maxr)
  
  if(method=="ed") {wsth=1
  }else{ wsth=max(1, sth-1)}
  
  for (r in wsth:sth){
    for (i in 1:ns){
      work= 1- xi*(xdat[i,r]-mut[i])/sigt[i]
      if(work < 1e-10 ) work=1e-10
      ztilda[i,r]= -log(work)/xi
    }
  }
  ztilda
}     
#------------------------------------
#-------------------------
ed.gum.gof=function(data,r=NULL,parx=NULL){
  
  R=r
  n=nrow(data)
  theta=parx   #theta[3]=-1e-4
  
  Diff <- ( eva::dgevr(data[, 1:R], loc = theta[1], scale = theta[2], 
                  shape = theta[3],log.d = TRUE)
            - eva::dgevr(data[, 1:(R -1)], loc = theta[1], 
                  scale = theta[2], shape = theta[3], 
                  log.d = TRUE) )
  
  EstVar <- sum((Diff - mean(Diff))^2)/(n - 1)
  FirstMom <- -log(theta[2]) -1 +(1 + theta[3])*digamma(R)
  Diff <- sum(Diff)/n
  Diff <- sqrt(n) * (Diff - FirstMom)/sqrt(EstVar)
  p.value <- 2 * (1 - pnorm(abs(Diff)))
  
  out <- list(statistics = as.numeric(Diff), 
              p.value = as.numeric(p.value), 
              theta = theta, 
              ybar = as.numeric(FirstMom))
  out
}
#---------------------------------
#------------------------------------------------------------
rcond.gof.gev11 = function(xdat=NULL, r=NULL, model='rgev11',
                           parx=NULL){
  
  ns= nrow(xdat)
  dim = ncol(xdat)
  
  y= rep(NA,ns)  
  ti= seq(1,ns)

  mu= parx[1]+ parx[2]*ti
  sig= exp(parx[3]+ parx[4]*ti)
  xi=parx[5]
  par11= cbind(mu, sig, rep(xi,ns))
  
  if(model=='rgev11'){ npar=5; mod ='gev'
  }else if(model=='rk4d11'){npar=6; mod='kap'}
  
  cdf.name=paste('cdf',mod,sep="")
  
  # Bartlett transformation to Uniform random variables
  
  sth=r
  for (t in 1:ns){
    
   if(sth==1) {
    
    y[t]= match.fun(cdf.name)(xdat[t,1], 
                              vec2par(par11[t,], mod,paracheck=F))
    
   }else if(sth >= 2){
    
    num= match.fun(cdf.name)(xdat[t,sth], 
                             vec2par(par11[t,], mod,paracheck=F))
    den= match.fun(cdf.name)(xdat[t,sth-1], 
                             vec2par(par11[t,], mod,paracheck=F))
    
    num[num==0]=1e-50
    den[den==0]=1e-50
    if(model=='rgev11'){
      a_s= 1
    }else if(model=='rk4d11'){
      a_s= 1-(sth-1)*parx[5]  
    }
    
    y[t]= exp(a_s*(log(num)-log(den))) 
    
   } # end if sth==1
    
  } #end for t
  
  sort(y)
}

#-------------------------------------------------------
spacing.gof.gev11 = function(xdat=NULL, r=NULL, parx=NULL){   
  
  ns= nrow(xdat); dim = ncol(xdat)
  D= rep(NA,ns)  
  ti= seq(1,ns)
  
  mu= parx[1]+ parx[2]*ti
  sig= exp(parx[3]+ parx[4]*ti)
  xi=parx[5]
  par11= cbind(mu, sig, rep(xi,ns))
  
  sth=r   # Tawn approach
  if(sth >= 1 & sth <= (dim-1)){  
      num = par11[,2] -par11[,3]*(xdat[,sth+1]-par11[,1])
      den = par11[,2] -par11[,3]*(xdat[,sth]-par11[,1])
      num[num <= 1e-50]=1e-50
      den[den <= 1e-50]=1e-50
      D= log(num/den)/par11[,3]
  } # end if sth
  sort(sth*D)
}
#-----------------------------------------------------

#------------------------------------
#-------------------------------------------------------------------
# ns=nrow(xdat)
# ff= (seq(1:ns)-0.35)/ns
# dim= ncol(xdat)
# 
# par(new=F); par(mfrow=c(2,5))
# lme.s=list()
# 
# for (sth in 1:dim){                        # sth=1
#   lme.s[[sth]] = parrgev.s(xdat,s=sth)
#   
#   if(lme.s[[sth]]$Ifail==0){
#     qua= quargev.s(ff, s=sth, para=lme.s[[sth]]$para)
#     
#     plot(sort(xdat[,sth]),qua,
#          main=paste("marg qq, s=",sth,sep=" "), pch=16)
#     abline(a=0,b=1)
#   }else{
#     cat("non-zero Ifail, s=",sth,"\n")
#   }
# }
# 
# cvm.test.lmomco(qua,lme.s[[dim]])