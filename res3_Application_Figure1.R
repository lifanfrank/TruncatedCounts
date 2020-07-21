#############################################
# source code to reproduce Figure 1
#############################################

setwd(".../ReproducibleResearch")  # need to change the directory here
source("calc.R")

###########################################
# a function to compute sample size 
# the estimated sample size is not restricted to an even number
###########################################

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
  ss=n0
  
  return(ss)
}    

### design parameters
Tset = c(Inf, seq(4,1,by=-1))   # truncation points
cvset=c(0,0.3,0.6,0.9)          # cv of cluster sizes
Nset_ind=Nset_asex=NULL

mbar=30                         # average cluster size
b0=log(2.7)+log(4/12)
b1=log(0.7)
s0=0.1                          # variance components in control arm
s1=0.1                          # variance components in intervention arm

# estimate sample size for sets of scenarios
set.seed(0611)
for(cv in cvset){
  temp1=temp2=NULL
  for(t in Tset){
    ccc<-calc(seed=0611,b0=b0,b1=b1,s0=s0,s1=s1,t=t)
    temp1<-c(temp1,sample_size(ccc[3],ccc[4],ccc[5],ccc[6],ccc[7],ccc[8],ccc[9],ccc[10],
                               mbar=mbar,cv=cv,cor="IND"))
    temp2<-c(temp2,sample_size(ccc[3],ccc[4],ccc[5],ccc[6],ccc[7],ccc[8],ccc[9],ccc[10],
                               mbar=mbar,cv=cv,cor="EXH"))
    print(t)
  }
  Nset_ind<-rbind(Nset_ind,temp1)
  Nset_asex<-rbind(Nset_asex,temp2)
}

# > Nset_ind
# [,1] [,2] [,3] [,4] [,5]
# temp1   39   39   40   44   62
# temp1   41   41   42   47   64
# temp1   48   48   49   53   71
# temp1   60   60   60   64   82
# > Nset_asex
# [,1] [,2] [,3] [,4] [,5]
# temp2   39   39   40   44   62
# temp2   40   40   41   45   63
# temp2   42   42   43   48   68
# temp2   46   47   48   54   77

library(plotrix)

# pdf(".../Fig1_application.pdf", width = 12, height = 6, paper = "special")  # change directory here
# par(mar=c(5.1, 6.1, 4.1, 2.1),mfrow=c(1,2))

###############################################################
# Panel (A) 
# corresponds to GEE analysis assuming independence working correlation
###############################################################

plot(c(1,3,4,5,6),Nset_ind[1,],col='blue',pch=16,ylim=c(39,82),
     xlab=expression(T),ylab=expression(hat(n)),xaxt="n",type='b',
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex=1.5,lwd=1.5,yaxt = 'n',
     main=expression(paste("(A)")))
axis(2,seq(40,80,by=5),las=1,cex.axis=1.5)
points(c(1,3,4,5,6),Nset_ind[2,],col='red',pch=17,type='b',cex=1.5,lwd=1.5)
points(c(1,3,4,5,6),Nset_ind[3,],col='brown',pch=18,type='b',cex=1.5,lwd=1.5)
points(c(1,3,4,5,6),Nset_ind[4,],col='black',pch=15,type='b',cex=1.5,lwd=1.5)

axis(1,at=c(1,3,4,5,6),labels=c(expression(infinity),4,3,2,1),
     cex.lab=1.5,cex.axis=1.5,)
axis.break(1,2,style="slash") 

legend(x="topleft",
       legend = c("CV = 0", "CV = 0.3", "CV = 0.6", "CV = 0.9"), 
       col=c("blue","red","brown","black"),lwd=2,cex=1.2,pch = c(16,17,18,15),
       bty='n')

###############################################################
# Panel (B) 
# corresponds to GEE analysis assuming arm-specific exchangeable working correlation
###############################################################

plot(c(1,3,4,5,6),Nset_asex[1,],col='blue',pch=16,ylim=c(39,82),
     xlab=expression(T),ylab=expression(hat(n)),xaxt="n",type='b',
     cex.lab=1.5,cex.axis=1.5,cex.main=1.5,cex=1.5,lwd=1.5,yaxt = 'n',
     main=expression(paste("(B)")))
axis(2,seq(40,80,by=5),las=1,cex.axis=1.5)
points(c(1,3,4,5,6),Nset_asex[2,],col='red',pch=17,type='b',cex=1.5,lwd=1.5)
points(c(1,3,4,5,6),Nset_asex[3,],col='brown',pch=18,type='b',cex=1.5,lwd=1.5)
points(c(1,3,4,5,6),Nset_asex[4,],col='black',pch=15,type='b',cex=1.5,lwd=1.5)

axis(1,at=c(1,3,4,5,6),labels=c(expression(infinity),4,3,2,1),
     cex.lab=1.5,cex.axis=1.5,)
axis.break(1,2,style="slash")

legend(x="topleft",
       legend = c("CV = 0", "CV = 0.3", "CV = 0.6", "CV = 0.9"), 
       col=c("blue","red","brown","black"),lwd=2,cex=1.2,pch = c(16,17,18,15),
       bty='n')

# dev.off()


