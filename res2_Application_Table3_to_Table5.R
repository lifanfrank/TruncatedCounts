###################################### 
# Source code to reproduce Table 3-5
######################################

setwd(".../ReproducibleResearch")  # need to change the directory here

library(xtable)
source("calc.R")

# function to explore right truncation on power
f.power<-function(mu0, mu1, tau0, tau1, kappa0, kappa1,rho0, rho1, n, mbar, cv, cor){
  
  # design parameters
  mu0<-as.numeric(mu0)
  mu1<-as.numeric(mu1)
  tau0<-as.numeric(tau0)
  tau1<-as.numeric(tau1)
  kappa0<-as.numeric(kappa0)
  kappa1<-as.numeric(kappa1)
  rho0<-as.numeric(rho0)
  rho1<-as.numeric(rho1)
  mbar<-as.numeric(mbar)
  cv<-as.numeric(cv)
  ss=n
  
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
  
  power=pt((qt(e1/2,df=ss-2) + abs(delta)/sqrt(sigma_b/ss)),df=ss-2)
  
  return(power)
}    

#################################################
# Table 3

# the power is the same under two working correlations
# under equal cluster sizes (cluster size cv = 0)
#################################################

###############################################
# First row
# # (e^{beta_0}, e^{\beta_1}) = (1.25,0.7)
###############################################
b0=log(1.25)
b1=log(0.7)
s2set=c(0.05,0.1,0.2,0.3,0.4)
nset=c(30,30,60,90,110)
mset=c(15,45,35,40,40)
Tset=c(Inf,6,5,4,3,2,1)
cv=0

FPOWER.ind=FPOWER.exh=matrix(NA,5,7)

set.seed(0611)
for(i in 1:5){
  for(j in 1:7){
    ccc<-calc(seed=0619,b0=b0,b1=b1,s0=s2set[i],s1=s2set[i],t=Tset[j])
    FPOWER.ind[i,j]<-f.power(ccc[3], ccc[4], ccc[5], ccc[6], ccc[7], ccc[8],ccc[9], ccc[10], n=nset[i], mbar=mset[i], cv=cv, cor="IND")
    FPOWER.exh[i,j]<-f.power(ccc[3], ccc[4], ccc[5], ccc[6], ccc[7], ccc[8],ccc[9], ccc[10], n=nset[i], mbar=mset[i], cv=cv, cor="EXH")
    print(j)
  }
}

# Latex code to generate Table entries
tab1=cbind(nset,mset,FPOWER.ind*100)
tab2=cbind(nset,mset,FPOWER.exh*100)  # no different from tab1
xtable(tab1,digits=1)

###############################################
# Second row
# # (e^{beta_0}, e^{\beta_1}) = (2.70,0.7)
###############################################
b0=log(2.70)
b1=log(0.7)
s2set=c(0.05,0.1,0.2,0.3,0.4)
nset=c(25,30,55,80,110)
mset=c(10,30,20,30,25)
Tset=c(Inf,6,5,4,3,2,1)
cv=0

FPOWER.ind=FPOWER.exh=matrix(NA,5,7)

set.seed(0611)
for(i in 1:5){
  for(j in 1:7){
    ccc<-calc(seed=0619,b0=b0,b1=b1,s0=s2set[i],s1=s2set[i],t=Tset[j])
    FPOWER.ind[i,j]<-f.power(ccc[3], ccc[4], ccc[5], ccc[6], ccc[7], ccc[8],ccc[9], ccc[10], n=nset[i], mbar=mset[i], cv=cv, cor="IND")
    FPOWER.exh[i,j]<-f.power(ccc[3], ccc[4], ccc[5], ccc[6], ccc[7], ccc[8],ccc[9], ccc[10], n=nset[i], mbar=mset[i], cv=cv, cor="EXH")
    print(j)
  }
}

# Latex code to generate Table entries
tab1=cbind(nset,mset,FPOWER.ind*100)
tab2=cbind(nset,mset,FPOWER.exh*100)    # no different from tab1
xtable(tab1,digits=1)

#################################################
# Table 4 and Table 5

# Unequal cluster sizes (cluster size cv = 0.6)
#################################################


###############################################
# First row
# # (e^{beta_0}, e^{\beta_1}) = (1.25,0.7)
###############################################
b0=log(1.25)
b1=log(0.7)
s2set=c(0.05,0.1,0.2,0.3,0.4)
nset=c(30,30,60,90,110)
mset=c(15,45,35,40,40)
Tset=c(Inf,6,5,4,3,2,1)
cv=0.6

FPOWER.ind=FPOWER.exh=matrix(NA,5,7)

set.seed(0611)
for(i in 1:5){
  for(j in 1:7){
    ccc<-calc(seed=0619,b0=b0,b1=b1,s0=s2set[i],s1=s2set[i],t=Tset[j])
    FPOWER.ind[i,j]<-f.power(ccc[3], ccc[4], ccc[5], ccc[6], ccc[7], ccc[8],ccc[9], ccc[10], n=nset[i], mbar=mset[i], cv=cv, cor="IND")
    FPOWER.exh[i,j]<-f.power(ccc[3], ccc[4], ccc[5], ccc[6], ccc[7], ccc[8],ccc[9], ccc[10], n=nset[i], mbar=mset[i], cv=cv, cor="EXH")
    print(j)
  }
}

# Latex code to generate Table entries
tab1=cbind(nset,mset,FPOWER.ind*100)   # Table 4: independence working correlation
tab2=cbind(nset,mset,FPOWER.exh*100)   # Table 5: arm-specific exchangeable working correlation
xtable(tab1,digits=1)
xtable(tab2,digits=1)

###############################################
# Second row
# # (e^{beta_0}, e^{\beta_1}) = (2.70,0.7)
###############################################
b0=log(2.70)
b1=log(0.7)
s2set=c(0.05,0.1,0.2,0.3,0.4)
nset=c(25,30,55,80,110)
mset=c(10,30,20,30,25)
Tset=c(Inf,6,5,4,3,2,1)
cv=0.6

FPOWER.ind=FPOWER.exh=matrix(NA,5,7)

set.seed(0611)
for(i in 1:5){
  for(j in 1:7){
    ccc<-calc(seed=0619,b0=b0,b1=b1,s0=s2set[i],s1=s2set[i],t=Tset[j])
    FPOWER.ind[i,j]<-f.power(ccc[3], ccc[4], ccc[5], ccc[6], ccc[7], ccc[8],ccc[9], ccc[10], n=nset[i], mbar=mset[i], cv=cv, cor="IND")
    FPOWER.exh[i,j]<-f.power(ccc[3], ccc[4], ccc[5], ccc[6], ccc[7], ccc[8],ccc[9], ccc[10], n=nset[i], mbar=mset[i], cv=cv, cor="EXH")
    print(j)
  }
}

# Latex code to generate Table entries
tab1=cbind(nset,mset,FPOWER.ind*100)      # Table 4: independence working correlation
tab2=cbind(nset,mset,FPOWER.exh*100)      # Table 5: arm-specific exchangeable working correlation
xtable(tab1,digits=1)
xtable(tab2,digits=1)



