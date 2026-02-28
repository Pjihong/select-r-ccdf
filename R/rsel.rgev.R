
# method="spacing"
# pplot=TRUE  # test="ad" 
# rsel.rgev(xdat, method="ed")
# rsel.rgev(xdat, method="ccdf")
# rsel.rgev(xdat, method="spacing")
#--------------------------------------------------------------------
rsel.rgev = function(xdat, sigL=0.05, num_inits=5, 
                     method=c("ed", "ccdf","spacing"),
                     test="cvm", pplot=FALSE){

  z=list();                                      
  ns=nrow(xdat); dim=ncol(xdat) 
  iD= mtheoU= matrix(NA, ns,dim)
  result = matrix(NA,dim,6)
  test.fun=paste(test,"test",sep=".")
  
  ED.result= gevrSeqTests.park(xdat, method='ed')
  
  if(method=="ed"){
    R=dim
    opt_r = min( which(ED.result$p.values <= sigL),R)
    if(opt_r==R){
      opt_r = min( which(ED.result$StrongStop <= sigL),R)
      if(opt_r==R){
       opt_r = min( which(ED.result$ForwardStop <= sigL),R)
      }
    }
  r.sel=opt_r
  } #end if method="ed"
  
  if(method != "ed"){

   sth=0
   while (sth < dim){      
     sth= sth+1 

     if(sth==1) {st.para=NULL
     }else{st.para= result[sth-1,2:4]}
     result[sth,1]= sth
    
     if(sth==1){ 
      gof = rgev.fit.park(xdat, r= sth, num_inits=num_inits, 
                         reltol=1e-5, start.para=st.para)
     
      if(gof$nllh >= 10^6){
       gof = rgev.fit.park(xdat, r= sth, num_inits=num_inits*2, 
                           reltol=1e-6, start.para=st.para)
      }
      result[sth,2:4]= gof$mle
      result[sth,5]= gof$nllh
      
     }else if(sth >= 2){
       
       if(ED.result$est.shape[sth-1] >= -0.5){
         
         result[sth,2]= ED.result$est.loc[sth-1]
         result[sth,3]= ED.result$est.scale[sth-1]
         result[sth,4]= ED.result$est.shape[sth-1]
         
       }else if(ED.result$est.shape[sth-1] < -0.5){ 
         
         gof = rk3d.fit.park(xdat, r= sth, h.fix= -0.001,
                             num_inits=num_inits, penk="CD",
                             reltol=1e-5)
         
         if(gof$nllh >= 10^6){
           gof = rgev.fit.park(xdat, r= sth, num_inits=num_inits*2, 
                               reltol=1e-6, start.para=st.para)
         }
         result[sth,2:4]= gof$mle[1:3]
         result[sth,5]= gof$nllh
       } # end if ED.result
     } # end if sth
     
     # cat("sth, mle=", sth, result[sth,2:4],"\n" )
     r.sel= dim+10
     
     if(method=="ccdf"){
       mtheoU[,sth] = rccdf.gof(xdat, r=sth, parx=result[sth,2:4]) 
       result[sth,6]= match.fun(test.fun)(mtheoU[,sth], null="punif")$p.value
      
#       cat("sth, p.val=", sth, result[sth,6],"\n" )
       
       if(result[sth,6] < sigL){r.sel= sth-1; break}
     }
     if(method=="spacing" & sth >= 2){
       iD[,sth]= spacing.gof(xdat, r= sth-1, parx=result[sth,2:4]) 
       result[sth,6]= match.fun(test.fun)(iD[,sth], null="pexp")$p.value 
       
#       cat("sth, p.val=", sth, result[sth,6],"\n" )
       
       if(result[sth,6] < sigL){r.sel= sth-1; break}
     }

  } #end while
} #end if method != "ed"
  
  if(r.sel== dim+10 & sth==dim) r.sel= dim
   work.r=min(r.sel+1,dim)

   z$method= method
   z$rsel= r.sel
   
   if( method != "ed"){
     colnames(result) = c("r","est.loc","est.scale","est.shape",
                        "nllh", "p.value")  #paste("p",method,sep="."))
     z$result = as.data.frame(result[1:work.r,])
     
   }else if(method=="ed"){
     z$result=ED.result
   }

#   png("D:pplot_spacing_101.png", width = 2400, height =2200, res = 350)
   
   if(pplot==TRUE){
     ff=(seq(1:ns)-0.35)/ns
     par(new=FALSE); par(mfrow=c(2,4))
     
     if(method=="ccdf"){
       for (sth in 1:work.r){
         plot(ff, mtheoU[,sth], main=paste("r=",sth,sep=" "),
              ylab="ccdf.model", 
              xlab=paste("EDF: p=",round(result[sth,6],3),sep=" "),pch=16) 
         abline(a=0,b=1) }
     }
     if(method=="spacing"){
       for(sth in 2:work.r){
         plot(ff, pexp(iD[,sth]), main=paste("r=",sth,sep=" "),
              ylab="spacing.model", 
              xlab=paste("EDF: p=",round(result[sth,6],3),sep=" "),
              pch=16, col="blue") 
         abline(a=0,b=1) }
     }
   } # end if pplot
   
#   dev.off()
 z
}
#-------------------------------------------------------
spacing.gof = function(xdat=NULL, r=NULL, parx=NULL){   
  
  ns= nrow(xdat); dim = ncol(xdat)

  sth=r                                         # Tawn approach
  if(sth >= 1 & sth <= (dim-1)){  # sth=9
      num = parx[2] -parx[3]*(xdat[,sth+1]-parx[1])
      den = parx[2] -parx[3]*(xdat[,sth]-parx[1])
      num[num <= 1e-50]=1e-50
      den[den <= 1e-50]=1e-50
      D= log(num/den)/parx[3]
  } # end if sth
 sort(sth*D)
}
#-------------------------------------------------------
#------------------------------------------------------------
rccdf.gof= function(xdat=NULL, r=NULL, parx=NULL){          
  
  ns= nrow(xdat); dim = ncol(xdat)
  cond.y= rep(NA,ns)  
  
  sth=r
  if(sth==1) {   
    cond.y[1:ns]= cdfgev(xdat[1:ns,1],    
                         vec2par(parx,'gev',paracheck=F))
      
  }else if(sth >= 2){ 
    num= cdfgev(xdat[,sth],              
                vec2par(parx,'gev',paracheck=F))
    den= cdfgev(xdat[,sth-1], 
                vec2par(parx,'gev',paracheck=F))
    num[num <= 1e-50]=1e-50
    den[den <= 1e-50]=1e-50

    cond.y[1:ns]= exp(log(num)-log(den))
  }
  sort(cond.y)
}
#----------------------------------------------------------
cdfrgev.s= function(x, s, para){
  w= (1-para[3]*(x-para[1])/para[2])^(1/para[3])
  1- pgamma(w, shape=s)
}
#-------------------------------------------------------------------
# ns=nrow(xdat)
# ff= (seq(1:ns)-0.35)/ns
# dim= ncol(xdat)
# 
# par(new=F); par(mfrow=c(2,5))
# lme.s=list()
# 
# for (sth in 1:dim){                        # sth=1
#   lme.s[[sth]] = parrgev.s(xdat,s=sth)     # lme for marginal dist
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
# cvm.test.lmomco(qua, lme.s[[dim]])