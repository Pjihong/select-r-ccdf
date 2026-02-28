# library(lmomco)
# 
# # Setting parameters for r-largest_GEV11 --------------------------------------
# 
# mu0= 100;
# mu1= 0.2
# sig0= log(10); 
# sig1= 0.02
# xi= -0.4  
# paraT=c(mu0, mu1, sig0, sig1, xi)
# # for real work, these parameters should be changed  
# 
# # Setting options for r-largest_GEV11 -----------------------------------------
# 
# nsample <- 60  # sample size (t)
# rlarg   <- 10   # r-largest
# bsample <- 1 # number of bootstrap
# 
# # Generate random number from r-largest_GEV11 ---------------------------------
# 
# rgev11_sample <-ran.rgev11(paraT,nsample,
#                            rlarg,bsample) # list format
# 
# xdat= rgev11_sample

#--------- Function for r-largest_GEV11 ------------------------------------------------
ran.rgev11 = function(para, nsample=30, rlarg=5, 
                      bsample=100){

  mu0 =para[1]
  mu1 =para[2]
  sig0=para[3]
  sig1=para[4]
  xi  =para[5]
  
  rand_array  <-array(rep(NA,nsample *rlarg *bsample),c(nsample, rlarg, bsample))
  rand_list   <-lapply(seq(bsample), function(x) rand_array[ , , x]) 
  
  umat_array  <-array(runif(nsample *rlarg *bsample,max=1,min=0),c(nsample, rlarg, bsample))
  umat_list   <-lapply(seq(bsample), function(x) umat_array[ , , x]) # array to list
  umat_rlag   <-lapply(seq(bsample), function(x) t(apply(umat_list[[x]],1,cumprod)))   # r-largest
  
  mut <- mu0+mu1*c(1:nsample)
  sigt<- exp(sig0+sig1*c(1:nsample))
  xi  <- xi
  
  savens <-lapply(seq(nsample), function(x) vec2par(c(mut[x],sigt[x],xi),'gev'))
    
  for(k in 1:bsample){
    for(t in 1:nsample){
      rand_list[[k]][t,]<-par2qua(umat_rlag[[k]][t,], savens[[t]], paracheck = FALSE)
    }
  }  
  return(rand_list)
}
#-----------------------------------------------------------