#------------------------------------------------------------
rccdf.gof.new= function(xdat=NULL, r=NULL, parx=NULL,
                    model="rgev"){          
  
  ns= nrow(xdat); dim = ncol(xdat)
  cond.y= rep(NA,ns)  
  
  sth=r
  if(model=="rgev"){
    
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
    
  }else if(model=="rk4d"){
    
    if(sth==1) {  
      cond.y[1:ns]= cdfkap(xdat[1:ns,1],    
                           vec2par(parx,'kap',paracheck=F))
      
    }else if(sth >= 2){ 
      num= cdfkap(xdat[,sth],              
                  vec2par(parx,'kap',paracheck=F))
      den= cdfkap(xdat[,sth-1], 
                  vec2par(parx,'kap',paracheck=F))
      num[num <= 1e-50]=1e-50
      den[den <= 1e-50]=1e-50
      
      cond.y[1:ns]= exp( (1-(sth-1)*parx[4])* (log(num)-log(den)) )
    }
  }
  sort(cond.y)
}
#----------------------------------------------------------
cdfrgev.s= function(x, s, para){
  w= (1-para[3]*(x-para[1])/para[2])^(1/para[3])
  1- pgamma(w, shape=s)
}
#-------------------------------------------------------------------