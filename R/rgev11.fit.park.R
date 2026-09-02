
# rgev11.fit.park(xdat,r=4,num_inits=10,lowb=-0.7,penk="CD")

# #--------------------------------------------------------------
# #-------------------------------------------------------    
# rgev11.fit.park = function(xdat=NULL, r=NULL, num_inits=15, 
#                          reltol=1e-6, show=FALSE, start.para=NULL,
#                          penk=NULL){
#   
#   rfit=list()  
#   dtr= xdat
#   ntry=num_inits
#   if(is.null(r)) r=ncol(xdat)
#   numr=r                              # r=numr=1
#   
#   # if(numr==1){
#   #   rfit= GEV110.park(data= as.vector(dtr[,1]), 
#   #                     ntry=num_inits, pretheta=start.para) 
#   # 
#   # }else{
#     rfit= rgev11.mle.park(xdat=dtr, numr=numr, 
#                           ntry=num_inits, lowb= -0.7, 
#                           reltol=reltol, const=TRUE,
#                           start.para=start.para,
#                           penk=penk)
#   # }
#   invisible(rfit)
#   return(rfit) 
# }

# T = "TRUE" const=T
# penk= "CD"

 # newdat=as.matrix(xdat)
 # rgev11.fit.park(newdat, r=1)
#------------------------------------------------------
#-------------------------------------------------------------  
rgev11.fit.park = function(xdat, r=NULL, num_inits=10, lowb= -1.0, 
                        reltol=1e-6, const=TRUE, maxit=200,
                        start.para=NULL, penk=NULL){
  
  k=list(); z=list()  # numr=4, ntry=6
 
  if(is.null(r)) r=ncol(xdat)    # r=numr=4 ;
  ntry=num_inits
  ns=nrow(xdat)
  xdatu <- xdat[, 1:r, drop = FALSE]

  ti <- seq(1:ns)       
  
  npar=5
  
  init= matrix(0, nrow=ntry, ncol=npar)
  init <- ginit.max.lme(data= xdat[,1], ntry, 
                        pretheta=start.para) 
  
  if(!is.null(start.para)) init[ntry,]=start.para
  nllh= rep(NA, ntry)
  
  tryCatch( 
    for(i in 1:nrow(init)){  # i=1
      
      work <-  try( optim(par=init[i,], fn=rlarg.lik.gev11, hessian = F, 
                          method=c("BFGS"), #method = c("Nelder-Mead"), 
                     #lower= c(-Inf, -Inf, -1.0), upper=c(Inf, Inf, -lowb), 
                     # coles style para
                     control = list(maxit = maxit, reltol=reltol),
                     lowb=lowb, const=const,
                     ti=ti, xdatu=xdatu, penk=penk, r=r),
                silent=T)  
                             # coles style para
      
      if(is(work)[1]=="try-error"){
        k[[i]] <- list(work=10^6)
      }else{
        k[[i]] <- work
        nllh[i]= k[[i]]$value
      }
      
    } #for  
  ) #tryCatch
  
  selc_num = which.min(nllh)

  z$conv <- k[[selc_num]]$convergence
  z$nllh <- min(nllh, na.rm=T)
  z$mle <- k[[selc_num]]$par             # coles style parameter
  z$mle[npar] = - k[[selc_num]]$par[npar]       # Hosking style para

  z
}
#-------------------------------------------------------------------- 
rlarg.lik.gev11 <- function(a=NULL, lowb=lowb, const=const,
                            ti=ti, xdatu=xdatu, 
                            penk=penk, r=r) {
  
  # mu0 <- a[1] #mulink(drop(mumat %*% (a[1:npmu])))
  # mu1 <- a[2] #siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
  # sig0 = a[3]
  # sig1 = a[4]
  
  mu= a[1]+a[2]*ti 
  sc= exp(a[3]+ a[4]*ti)
  xi <- a[5] # coles style

  if(const==TRUE & xi < lowb) return(10^6)
  
  if (any(sc <= 0)) return(10^6)
  
  y <- 1 + xi * (xdatu - mu)/sc
  
  if (min(y, na.rm = TRUE) <= 0) 
    return(10^6)
  else {
    y <- (1/xi + 1) * log(y) + log(sc)
    y <- rowSums(y, na.rm = TRUE)
    
    l <- sum( (1 + xi * (xdatu[,r] - mu)/sc)^(-1/xi) + y)  # xdatu[,r]
  }
  
  if (is.null(penk) == FALSE) {
    if (penk == "CD") {
      if (xi[1] >= 0) {
        p_k <- 1
      }
      else if (xi[1] > -1 && xi[1] < 0) {
        p_k <- exp(-((1/(1 + xi[1])) - 1))
      }
      else if (xi[1] <= -1) {
        p_k <- 1e-50
      }
    }
    else if (penk == "MS") {
      if (xi[1] >= 0.5 || xi[1] <= -0.5) {p_k =1e-50
      }else{ p_k <- ((0.5+ xi[1])^(6-1)) * ((0.5- xi[1])^(9- 
                                                             1))/beta(6, 9)}
    }
  } #end if is.null
  
  if(is.null(penk)) {penalty =0
  }else{ penalty <- r * log(p_k) }
  
  l -penalty  
}
#-----------------------------------------------
#------------------------------------------------------
ginit.max.lme <-function(data,ntry=ntry, pretheta=NULL){
  
  init <-matrix(0, nrow=ntry, ncol=5)
  
  lmom_init = lmoms(data,nmom=5)
  lmom_est <- pargev(lmom_init)
  
  init[1,1] <- lmom_est$para[1]
  init[1,2] <- 0.001
  init[1,3] <- log(lmom_est$para[2])
  init[1,4] =  -0.001
  init[1,5] = lmom_est$para[3]
  
  if(is.null(pretheta)) {
    maxm1=ntry; maxm2=ntry-1
  }else{
    maxm1=ntry-1; maxm2=ntry-2
    init[ntry,1:5]= pretheta[1:5]}
  
  init[2:maxm1,1] <- init[1,1]+rnorm(n=maxm2,mean=0,sd = 20)
  init[2:maxm1,2] <- rnorm(n=maxm2,mean=0,sd = 1)
  init[2:maxm1,3] <- init[1,3]+rnorm(n=maxm2,mean=0,sd = 2)
  init[2:maxm1,4] <- rnorm(n=maxm2,mean=0,sd = 0.5)
  init[2:maxm1,5] <- runif(n=maxm2,min= -0.5,max=0.5)
  
  return(init)
}
#----------------------------------------------------
# 
# #---------------------------------------------    
# rgev.lik.cvnll.park <- function(a, r=NULL, xdat=NULL) {
#   mu <- a[1] # mulink(drop(mumat %*% (a[1:npmu])))
#   sc <- a[2] # siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
#   
#   xi <- -a[3] # hosking style xi         # shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))
#   
#   if (any(sc <= 0)) 
#     return(10^6)
#   
#   xdatu <- xdat[, 1:r, drop = FALSE]
#   u <- apply(xdatu, 1, min, na.rm = TRUE)
#   
#   y <- 1 + xi * (xdatu - mu)/sc            # xi: coles style
#   if (min(y, na.rm = TRUE) <= 0) {
#     l <- 10^6
#   }else {
#     y <- (1/xi + 1) * log(y) + log(sc)
#     y <- rowSums(y, na.rm = TRUE)
#     l <- sum((1 + xi * (u - mu)/sc)^(-1/xi) + y)
#   }
#   l
# }
# #-------------------------------------------------------   
# #--------------------------------------------------------------
# #  The following function is just a simple modification of 
# #  'rlarg.fit' function in the 'ismev' package.
# #--------------------------------------------------------------
# 
# rlarg.fit.gev11 = function (xdat, r =NULL, start=NULL, ydat = NULL, mul = NULL, sigl = NULL, 
#                         shl = NULL, mulink = identity, siglink = identity, shlink = identity, 
#                         muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, 
#                         method = "Nelder-Mead", maxit = 2000, lowb=-1, const=NULL, 
#                         reltol= 1e-6, penk=NULL, ...) 
# {
# 
#   z <- list()     
#   if(is.null(r)) r=ncol(xdat)    # r=numr=4 ; init=init[1,]
#   data1= xdat[,1]                                # coles style para
#   
#   npmu <- length(mul) + 1
#   npsc <- length(sigl) + 1
#   npsh <- length(shl) + 1
#   z$trans <- FALSE
#   in2 <- sqrt(6 * var(xdat[, 1]))/pi
#   in1 <- mean(xdat[, 1]) - 0.57722 * in2
#   
#   # if (is.null(mul)) {
#   #   mumat <- as.matrix(rep(1, dim(xdat)[1]))
#   #   if (is.null(muinit)) 
#   #     muinit <- in1
#   # }  else {
#   #   z$trans <- TRUE
#   #   mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
#   #   if (is.null(muinit)) 
#   #     muinit <- c(in1, rep(0, length(mul)))
#   # }
#   # if (is.null(sigl)) {
#   #   sigmat <- as.matrix(rep(1, dim(xdat)[1]))
#   #   if (is.null(siginit)) 
#   #     siginit <- in2
#   # }  else {
#   #   z$trans <- TRUE
#   #   sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
#   #   if (is.null(siginit)) 
#   #     siginit <- c(in2, rep(0, length(sigl)))
#   # }
#   # if (is.null(shl)) {
#   #   shmat <- as.matrix(rep(1, dim(xdat)[1]))
#   #   if (is.null(shinit)) 
#   #     shinit <- 0.1
#   # }  else {
#   #   z$trans <- TRUE
#   #   shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
#   #   if (is.null(shinit)) 
#   #     shinit <- c(0.1, rep(0, length(shl)))
#   # }
#   xdatu <- xdat[, 1:r, drop = FALSE]
#   
#   # if(is.null(init)) {init <- c(muinit, siginit, shinit)    # park modified these 2 lines
#   # }else{init = init}
#   init=start
#   
#   z$model <- list(mul, sigl, shl)
#   z$link <- deparse(substitute(c(mulink, siglink, shlink)))
#   u <- apply(xdatu, 1, min, na.rm = TRUE)
#   
#   ti <- matrix(nrow=length(data1),ncol=1)
#   ti[,1] <- seq(1:length(data1))
#   
#   # if(is.null(const) | const==FALSE){
#   #   
#   #     x <- optim(init, rlarg.lik, hessian = F, method = method, 
#   #              control = list(maxit = maxit, reltol=reltol,...))
#   # }else{
#   # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
# 
#      x <- optim(init, rlarg.lik.gev11, hessian = F, #method = c("L-BFGS-B"), 
#              #lower= c(-Inf, -Inf, -1.0), upper=c(Inf, Inf, -lowb),     # coles style para
#              control = list(maxit = maxit, reltol=reltol),
#              ti=ti, xdatu=xdatu, u=u, penk=penk)
#   # }
#   #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
#   #+
#   # mu0 <- x$par[1] #mulink(drop(mumat %*% (a[1:npmu])))
#   # mu1 <- x$par[2] #siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
#   # sig0 = x$par[3]
#   # sig1 = x$par[4]
#   # xi <- x$par[5]   
#   # mu <- mulink(drop(mumat %*% (x$par[1:npmu])))
#   # sc <- siglink(drop(sigmat %*% (x$par[seq(npmu + 1, length = npsc)])))
#   # xi <- shlink(drop(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)])))
#   
#   
#   z$conv <- x$convergence
#   z$nllh <- x$value
#   # #  z$data <- xdat
#   # if (z$trans) {
#   #   for (i in 1:r) z$data[, i] <- -log((1 + (as.vector(xi) * 
#   #                                              (xdat[, i] - as.vector(mu)))/as.vector(sc))^(-1/as.vector(xi)))
#   # }
#   z$mle <- x$par                       # coles style parameter
#   # z$cov <- solve(x$hessian)
#   # z$se <- sqrt(diag(z$cov))
#   #  z$vals <- cbind(mu, sc, xi)
#   z$r <- r
#   # if (show) {
#   #   if (z$trans) 
#   #     print(z[c(2, 3)])
#   #   print(z[4])
#   #   if (!z$conv) 
#   #     print(z[c(5, 7, 9)])
#   # }
#   class(z) <- "rlarg.fit.gev11"
#   invisible(z)
#   
#   return(z)
# }
# # ----------------------------------------------------------------
