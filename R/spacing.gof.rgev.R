#-------------------------------------------------------
spacing.gof = function(xdat=NULL, r=NULL, model='rgev',
                       parx=NULL){       
  
  ns= nrow(xdat)
  dim = ncol(xdat)
  D= rep(NA,ns)  
  
  sth=r                                         # Tawn approach
  if(sth >= 1 & sth <= (dim-1)){ 
    num = parx[2] -parx[3]*(xdat[,sth+1]-parx[1])
    den = parx[2] -parx[3]*(xdat[,sth]-parx[1])
    num[num==0]=1e-50
    den[den==0]=1e-50
    D= log(num/den)/parx[3]
  } # end if sth
  sort(sth*D)
}

#-------------------------------------------------------