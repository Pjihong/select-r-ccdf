library(eva)
library(goftest)
library(evmr)
library(ismev)
library(lmomco)
library(Rsolnp)

# xdat= na.omit(sancheong[,2:21])
# s1=Sys.time()
#  rsel.rgev.SeqStop(xdat, sigL=0.05, method="ed")
# Sys.time()-s1
# 
# s1=Sys.time()
#  rsel.rgev.SeqStop(xdat, sigL=0.05, method="ccdf")
# Sys.time()-s1
# 
# s1=Sys.time()
#   rsel.rgev.SeqStop(xdat, sigL=0.05, method="spacing", pplot=TRUE)
# Sys.time()-s1


#--------------------------------------------------------------------
rsel.rgev.SeqStop = function(xdat, sigL=0.05, num_inits=10, 
                           method=c("ed", "ccdf","spacing"),
                           test="cvm", pplot=FALSE){
  
  z=list();                                      
  ns=nrow(xdat); dim=ncol(xdat) 
  iD= mtheoU= matrix(NA, ns,dim)
  sp.result= result = matrix(NA,dim,8)
  ed1 = matrix(NA,1,8)
  test.fun=paste(test,"test",sep=".")
  
  ED0= gevrSeqTests.park(xdat, method='ed')
  
  colnames(ed1)= colnames(ED0)
  
  gev1= rgev.fit.park(xdat, r= 1, num_inits=num_inits, 
                      reltol=1e-5)

  ed1[1,1] =1
  Ugev = cdfgev(sort(xdat[,1]), 
                vec2par(gev1$mle,"gev") )
  
  testU01= match.fun(test.fun)(Ugev, null="punif")
  
  ed1[1,2]= as.numeric(testU01$p.value)  
  ed1[1,5]= as.numeric(testU01[1])      # stat
  ed1[1,6:8]= as.numeric(gev1$mle)      # para est
  
  ED.result= rbind(ed1, ED0)
  
  # ED.result[, 3] <- rev(pSeqStop(rev(ED.result[, 2]))$ForwardStop)
  # ED.result[, 4] <- rev(pSeqStop(rev(ED.result[, 2]))$StrongStop)
  
  if(method=="ed"){
    R=dim
    opt_r = min( which(ED.result$p.values <= sigL),R)
    if(opt_r==R){
      opt_r = min( which(ED.result$ForwardStop <= sigL),R)
      if(opt_r==R){
        opt_r = min( which(ED.result$StrongStop <= sigL),R)
      }
    }
    
    if(any(ED.result$p.values <= sigL) |
       any(ED.result$ForwardStop <= sigL) |
       any(ED.result$StrongStop <= sigL)) { r.sel= opt_r-1
    }else{ r.sel=R }
    
  } #end if method="ed"
  
  if(method != "ed"){
    
    for (sth in 1:dim){      
      
      if(sth==1) { st.para= NULL
      }else{       st.para= result[sth-1,6:8] }
      result[sth,1]= sth
      
      if(sth==1){ 

        result[1,6:8]= gev1$mle
        result[1,5]= gev1$nllh
        
      }else if(sth >= 2){
        
        if(ED.result$est.shape[sth] >= -0.5){
          
          result[sth,6]= ED.result$est.loc[sth]
          result[sth,7]= ED.result$est.scale[sth]
          result[sth,8]= ED.result$est.shape[sth]  # hosking style
          
        }else if(ED.result$est.shape[sth] < -0.5){ 
          
          gof = rk3d.fit.park(xdat, r= sth, h.fix= -0.001,
                              num_inits=num_inits, penk="CD",
                              reltol=1e-5)
          
          if(gof$nllh >= 10^6){
            gof = rgev.fit.park(xdat, r= sth, 
                                num_inits=num_inits*2, 
                                reltol=1e-6, 
                                start.para=st.para)
          }
          result[sth,6:8]= gof$mle[1:3]
          result[sth,5]= gof$nllh
        } # end if ED.result
        
      } # end if sth

      if(method=="ccdf"){
        
        mtheoU[,sth] = rccdf.gof.new(xdat, r=sth, 
                             parx=result[sth,6:8], model="rgev") 
        
        result[sth,2]= match.fun(test.fun)(mtheoU[,sth], 
                                           null="punif")$p.value
      }
      
      if( method=="spacing" & sth == 1){

        result[1,2]= match.fun(test.fun)(Ugev, 
                                  null="punif")$p.value
      }
      if( method=="spacing" & sth >= 2){
        
        iD[,sth]= spacing.gof(xdat, r= sth-1, 
                              parx=result[sth,6:8]) 
        
        result[sth,2]= match.fun(test.fun)(iD[,sth], 
                                    null="pexp")$p.value 
      }
      
    } #end for sth

    result[, 3] <- rev(pSeqStop(rev(result[, 2]))$ForwardStop)
    result[, 4] <- rev(pSeqStop(rev(result[, 2]))$StrongStop)
    
    colnames(result) =  colnames(ED0)

    result= as.data.frame(result)
    
      R=dim
      opt_r = min( which(result$p.values <= sigL),R)
      if(opt_r==R){
        opt_r = min( which(result$ForwardStop <= sigL),R)
        if(opt_r==R){
          opt_r = min( which(result$StrongStop <= sigL),R)
        }
      }
      
      if(any(result$p.values <= sigL) |
         any(result$ForwardStop <= sigL) |
         any(result$StrongStop <= sigL)) { r.sel= opt_r-1
      }else{ r.sel=R }

  } #end if method != "ed"
  
  work.r=dim
  
  z$method= method
  z$rsel= r.sel
  
  if( method != "ed"){

    z$result = as.data.frame(result[1:work.r,])
    
  }else if(method=="ed"){

    z$result= ED.result
  }

  if(pplot==TRUE){
    ff=(seq(1:ns)-0.35)/ns
    par(new=FALSE); par(mfrow=c(2,4))
    
    if(method=="ccdf"){
      for (sth in 1:work.r){
        plot(ff, mtheoU[,sth], main=paste("r=",sth,sep=" "),
             ylab="ccdf.model", 
             xlab=paste("EDF: p=",round(result[sth,2],3),sep=" "),
             pch=16) 
        abline(a=0,b=1) }
    }
    if(method=="spacing"){
      
      plot(ff, Ugev,
           main=paste("r=",1,sep=" "),
           ylab="spacing.model", 
           xlab=paste("EDF: p=",round(result[1,2],3),sep=" "),
           pch=16, col="blue") 
      abline(a=0,b=1)
      
      for(sth in 2:work.r){
        plot(ff, pexp(iD[,sth]), main=paste("r=",sth,sep=" "),
             ylab="spacing.model", 
             xlab=paste("EDF: p=",round(result[sth,2],3),sep=" "),
             pch=16, col="blue") 
        abline(a=0,b=1) }
    }
  } # end if pplot
  z
}
#-------------------------------------------------------