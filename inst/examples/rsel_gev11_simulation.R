library(selectrgev)


#------------------------------------------------------
truer=rcontam= 4;
maxr= 8; 
sigL=sigL_ed= 0.05
sigL_cvm= 0.05
ntest= 3               
model=rLOS= "rgev11"

maxrun= 200  # number of trials = 200

kpa=c(-0.35, -0.2, -1e-4, 0.2, 0.35)

nsam=c(30,50,80) # sample size =c(30, 50, 80)


for (insam in 1:length(nsam)){
  
  n=nsam[insam]

for (ipa in 1:length(kpa)){
  
  para=c(0, 0.1, 1, .02, kpa[ipa])  
  # parameter=(mu0, mu1, sig0, sig1, k)

r_est=matrix(NA, maxrun, ntest)
Bader=list()

for (irun in 1:maxrun){  
  
  Bader[[irun]]= contam_rlosr(n=n, maxr=maxr, model=rLOS,
                              para=para,
                              rcontam=rcontam, mixp=0.5)$cdata
  
  xdat=as.matrix(Bader[[irun]])
  
  # Edtest= multi.rEdtest.park(xdat, model=rLOS, method='ed')
  # r_est[irun,1]= one_optr(Edtest, mid=mid, sigL=sigL_ed, h.fix=NULL)
  
  r_est[irun,1]= rsel.rgev11(xdat, model='rgev11',
                             method="ed", seq.cut=TRUE,
                             qqplot=FALSE, sigL=sigL_cvm)$r.sel
  
  r_est[irun,2]= rsel.rgev11(xdat, model='rgev11',
                             method="ccdf", seq.cut=TRUE,
                             qqplot=FALSE, sigL=sigL_cvm)$r.sel
  
  r_est[irun,3]= rsel.rgev11(xdat, model='rgev11',
                             method="spacing", seq.cut=TRUE,
                             qqplot=FALSE, sigL=sigL_cvm)$r.sel
  
  cat("irun, truer, r.cvm=",
      irun, truer, r_est[irun,],"\n")
  
}  #end for irun


bias= rmse= se= rep(NA,ntest)

for (k in 1:ntest){
  bias[k]= mean(r_est[,k]- truer, na.rm=TRUE)
  rmse[k]= sqrt( bias[k]^2 + var(r_est[,k]) )
}
se= sqrt(rmse^2 - bias^2)
cat("bias,se,rmse=",bias, se, rmse,"\n")


} # end for ip
  
} #end for in

