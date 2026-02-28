#-------------------------------------------------------   
gevrSeqTests.park = function (data, bootnum = NULL, method = c("ed", "pbscore", "multscore"), 
                              information = c("expected", "observed"), allowParallel = FALSE, 
                              numCores = 1,
                              low.xi= -0.7) 
{
  data <- as.matrix(data)
  R <- ncol(data)
  method <- match.arg(method)
  
  if (method != "ed") {
    if (is.null(bootnum)) 
      stop("Must enter the number of bootstrap replicates!")
    information <- match.arg(information)
    result <- matrix(0, R, 8)
    for (i in 1:R) {
      result[i, 1] <- i
      if (method == "multscore") 
        fit <- eva::gevrMultScore(data[, 1:i], bootnum, information)
      if (method == "pbscore") 
        fit <- eva::gevrPbScore(data[, 1:i], bootnum, information, 
                           allowParallel, numCores)
      result[i, 2] <- fit$p.value
      result[i, 5] <- fit$statistic
      # result[i, 6:8] <- fit$theta
      result[i, 6:7] <- fit$theta[1:2]
      result[i, 8] <-  -fit$theta[3]   #----- hosking style
    }
  }
  else {
    if (R == 1) 
      stop("R must be at least two")
    result <- matrix(0, R - 1, 8)
    for (i in 2:R) {
      result[i - 1, 1] <- i
      fit <- gevrEd.park1(data[, 1:i])
      
      if(fit$ifail==1){
        cat("gevrEd.park1: ifail, theta=", fit$ifail, fit$theta,"\n")
        mle <- rk3d.fit.park(data[, 1:i], h.fix=0.001, penk="CD")$mle
        
        mle[3]= -mle[3]
        fit = gevrEd.park1(data[, 1:i], theta=mle[1:3])
      }
      result[i - 1, 2] <- fit$p.value
      result[i - 1, 5] <- fit$statistic
      result[i - 1, 6:7] <- fit$theta[1:2]
      result[i - 1, 8] <-  -fit$theta[3]   #----- hosking style
    }
  }
  result[, 3] <- rev(pSeqStop(rev(result[, 2]))$ForwardStop)
  result[, 4] <- rev(pSeqStop(rev(result[, 2]))$StrongStop)
  colnames(result) <- c("r", "p.values", "ForwardStop", "StrongStop", 
                        "statistic", "est.loc", "est.scale",
                        "est.shape")
  as.data.frame(result)
}
#--------------------------------------------------------
gevrEd.park1 =function (data, theta = NULL) 
{
  data <- as.matrix(data)
  R <- ncol(data)
  if (R == 1) 
    stop("R must be at least two")
  n <- nrow(data)
  ifail=0
  
  if (is.null(theta)) {
    y <- tryCatch(eva::gevrFit(data, method = "mle"), error = function(w) {
      return(NULL)
    }, warning = function(w) {
      return(NULL)
    })
    if (is.null(y)) {
      warning("Maximum likelihood failed to converge at initial step in gevrFit")
      out=list(); out$ifail= 1
      
      y= rgev.fit.park(data)  # park
      theta = y$mle           # park  y$mle is hosking style
      theta[3]= -theta[3]

    }else{
      theta <- y$par.ests  # coles style
      #      theta[3]= -theta[3]  # hosking style
    }
  }
  
  if(theta[2] < 1e-4) theta[2]= 1e-4
  Diff1 <- eva::dgevr(data[, 1:R], loc = theta[1], scale = theta[2], 
                 shape = theta[3], log.d = TRUE)
  
  Diff2= eva::dgevr(data[, 1:(R-1)], loc = theta[1], scale = theta[2], 
               shape = theta[3], log.d = TRUE)
  
  Diff= Diff1 - Diff2
  
  EstVar <- sum((Diff - mean(Diff))^2)/(n - 1)
  FirstMom <- -log(theta[2]) - 1 + (1 +theta[3]) * digamma(R)
  Diff <- sum(Diff)/n
  Diff <- sqrt(n) * (Diff - FirstMom)/sqrt(EstVar)
  p.value <- 2 * (1 - pnorm(abs(Diff)))
  
  out <- list(as.numeric(Diff), as.numeric(p.value), theta, ifail)
  names(out) <- c("statistic", "p.value", "theta", "ifail")
  return(out)
}
