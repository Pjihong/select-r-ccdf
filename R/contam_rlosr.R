#---------------------------------------------------------
contam_rlosr = function(n=NULL, maxr=NULL, model=NULL, para=NULL,
                             rcontam=NULL, mixp=0.7){
  
  # para is hosking style 
  # n=30; maxr=10; para=c(100,10, -.3, -.5)
  # rcontam= 5; model="rk4d"

  z=list()
  rdata= matrix(NA,n,maxr)
  if(rcontam > maxr) stop("rcontam should not be greater than maxr")
  
  if(model=="rgev"){
     rdata= eva::rgevr(n,r=maxr,loc=para[1],scale=para[2], 
                  shape= -para[3])
  }
  if(model=="rk4d"){
    rdata= rk4dr(n,r=maxr,loc=para[1],scale=para[2], 
                 shape= para[3], shape2=para[4])$rmat
  }
  if(model=="rwak"){
    rdata= rwakr(n,r=maxr,loc=para[1],scale=para[2], 
                 shape= para[3], shape2=para[4])$rmat
  }
  if(model=="rgev11"){
    rdata = ran.rgev11(para=para, nsample=n, rlarg=maxr)
    # para=c(mu0, mu1, sig0, sig1, xi) 
    # where sig=exp(sig0 +sig1*t),
    rdata= as.matrix(rdata[[1]])
  }
  
  cdata=rdata
  
  if(rcontam <= maxr-2){
    
   for (rc in rcontam:(maxr-2)){
    cdata[,rc+1]= rdata[,rc+1]*mixp + rdata[,rc+2]*(1-mixp) 
   }
  } # end if
  
  # if(rcontam==maxr-1){
  #   cdata[,maxr]= rdata[,maxr-1]*mixp + rdata[,maxr]*(1-mixp)
  # }
  
  z$rdata= rdata; z$cdata= cdata
  z$rtrue= rcontam; z$mixp= mixp
  return(z)
}
#----------------------------------------------------------
rk4dr = function (n, r, loc = 0, scale = 1, shape1 = 0.1, shape2 = 0.1) 
{
  # the same code from the library evmr3
  
  z <- list()
  umat <- matrix(0, nrow = n, ncol = r)
  wmat <- matrix(0, nrow = n, ncol = r)
  colnames(umat) <- paste0("u", 1:r)
  colnames(wmat) <- paste0("w", 1:r)
  i <- 1
  while (i <= n) {
    u <- stats::runif(r)
    w <- numeric(r)
    w[1] <- u[1]
    for (j in 2:r) {
      w[j] <- w[j - 1] * (u[j]^(1/(1 - (j - 1) * shape2)))
    }
    if (all(diff(w) < 0)) {
      umat[i, ] <- u
      wmat[i, ] <- w
      i <- i + 1
    }
  }
  z$umat <- umat
  z$wmat <- wmat
  z$rmat <- qk4d(wmat, loc = loc, scale = scale, shape1 = shape1, 
                 shape2 = shape2)
  colnames(z$rmat) <- paste0("r", 1:r)
  invisible(z)
}
#--------------------------------------------------------
qk4d = function (p, loc = 0, scale = 1, shape1 = 0.1, shape2 = 0.1) 
{
  # the same code from the library evmr3
  
  mu = loc
  sig = scale
  xi = shape1
  h = shape2
  yp = (1 - (p)^h)/h
  z <- mu + (sig/xi) - (sig/xi) * (yp^xi)
  z
}
#----------------------------------------------------------
rwakr <- function(gen=100, n, r = maxr, para) {
  
  # para = re-parameterized 5 paras (Hosking and Wallis) 
  # generate r-largest for a row
  wakeby_generate_max_values <- function(gen, max_count, par) {
    # Uniform -> Wakeby quantile
    data <- quawak(f = runif(gen),
                   wakpara = vec2par(par, type = "wak", paracheck = FALSE),
                   paracheck = FALSE)
    sort(data, decreasing = TRUE)[1:max_count]
  }
  
  # make a matrix of n rows and r columns
  out <- t(replicate(n, wakeby_generate_max_values(gen, r, para)))
  #colnames(out) <- paste0("r", seq_len(r))
  out
}