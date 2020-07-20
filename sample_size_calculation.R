##########################################################################################
# Sample size new with constant cluster size and varying cluster sizes

# Input parameters:
# mu0, mu1: marginal means for two arms
# tau0, tau1: marginal variances for two arms
# kappa0, kappa1: marginal coefficient of variation of the outcome for the two arms
#                 kappa0=sqrt(tau0)/mu0
#                 kappa0=sqrt(tau1)/mu1
# rho0, rho1: marginal ICC for the two arms
# cv: coefficient of variation of the cluster size distribution
# cor: correlation structure
#      "IND" = independence
#      "EXH" = arm-specific exchangeable 
##########################################################################################

sample_size<-function(mu0, mu1, tau0, tau1, kappa0, kappa1,rho0, rho1, mbar, cv, cor){
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
  
  return(ss)
}    
    


