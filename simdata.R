##################################################
# Simulate correlated count data from a 
# Poisson mixed model allowing for arm-specific
# variance components

# n: number of clusters
# mbar: mean cluster sizes
# cv: cv of cluster sizes
# b0: beta0 in Poisson mixed model
# b1: beta1 in Poisson mixed model
# s0: variance component in control arm
# s1: variance component in intervention arm
# t: truncation point
##################################################

simdata=function(n,mbar,cv,b0,b1,s0,s1,t){
  require(truncdist)
  
  # cluster sizes
  if(cv==0){
    m=rep(mbar,n)
  } else{
    m=round(rgamma(n,shape=1/(cv^2),rate=1/(mbar*cv^2)))
    m[m<2]=2
  }
  N=sum(m)
  
  # treatment
  trt=rep(0,n)
  id=sample(1:n,size=n/2,replace=FALSE,prob=rep(0.5,n))
  trt[id]=1
  x=rep(trt,m)
  
  # cluster id
  cluster=rep(1:n,m)
  
  # poisson mixed model
  u0=rep(rnorm(n,0,sqrt(s0)),m)
  u1=rep(rnorm(n,0,sqrt(s1)),m)
  lambda=exp(b0+b1*x+u0*(1-x)+u1*x)
  y=NULL
  for(k in 1:N){
    y=c(y,rtrunc(1,spec="pois",a=-Inf,b=t,lambda=lambda[k]))
  }
  
  # combine dataset
  datatest=data.frame(y=y,x=x,cluster=cluster)
  return(list(datatest=datatest,trt=trt))
}



