# Sample size new with constant cluster size and varying cluster sizes following eq1 and eq3

sample_size_tp<-function(mu0, mu1, tau0, tau1, kappa0, kappa1,rho0, rho1, mbar, cv, cor){
  #set parameters
  e1=0.05 #type-1 error
  e2=0.2 #type-2 error
  
  delta<-log(mu1/mu0)
  pi=1/2 

  if (cor=="IND"){
    sigma_b<-kappa0^2/(1-pi)/mbar*(1+rho0*((1+cv^2)*mbar-1))+kappa1^2/pi/mbar*(1+rho1*((1+cv^2)*mbar-1))
  }  
  
  if (cor=="EXH"){
    sigma_b<-kappa0^2*(1+(mbar-1)*rho0)/(1-pi)/mbar*(1-cv^2*mbar*rho0*(1-rho0)/(1+(mbar-1)*rho0)^2)^-1 + kappa1^2*(1+(mbar-1)*rho1)/(1-pi)/mbar*(1-cv^2*mbar*rho1*(1-rho1)/(1+(mbar-1)*rho1)^2)^-1
  }
  
  n0=3
  n=(qt(e1/2,df=n0-2)+qt(e2,df=n0-2))^2*sigma_b/delta^2 
  
  while(n0<n){
    n=(qt(e1/2,df=n0-2)+qt(e2,df=n0-2))^2*sigma_b/delta^2
    n0=n0+1
  }
  if (n0%%2==1) ss=n0+1 else ss=n0
  
  power=pt(as.numeric((qt(e1/2,df=ss-2) + abs(delta)/sqrt(sigma_b/ss))),df=ss-2)

  return(c(round(ss,0),power))
}    
    


