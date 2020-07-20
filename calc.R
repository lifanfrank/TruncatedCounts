############################################################
# calculate marginal values based on Poisson mixed models
# this program implements the set of conversion formulas
# to obtain the design parameters required by the sample size formulas

# seed: set seed for Monte Carlo integration (with respect to gaussian random effects)
# b0: beta0 in Poisson mixed model
# b1: beta1 in Poisson mixed model
# s0: variance component in control arm
# s1: variance component in intervention arm
# t: truncation point
##############################################################

calc=function(seed=NULL,b0,b1,s0,s1,t){
  Q=function(t,eta,u){
    if(t>=0){
      lambda=exp(eta+u)
      return(ppois(t,lambda)*exp(lambda))
    }else{
      return(0)
    }
  }
  
  # mean
  fmu0=function(u){
    q0=Q(t,b0,u)
    q1=Q(t-1,b0,u)
    lambda=exp(b0+u)
    (q1/q0)*lambda
  }
  fmu1=function(u){
    q0=Q(t,b0+b1,u)
    q1=Q(t-1,b0+b1,u)
    lambda=exp(b0+b1+u)
    (q1/q0)*lambda
  }
  
  # second moment
  fmm0=function(u){
    q0=Q(t,b0,u)
    q1=Q(t-1,b0,u)
    q2=Q(t-2,b0,u)
    lambda=exp(b0+u)
    (lambda^2*q2+lambda*q1)/q0
  }
  fmm1=function(u){
    q0=Q(t,b0+b1,u)
    q1=Q(t-1,b0+b1,u)
    q2=Q(t-2,b0+b1,u)
    lambda=exp(b0+b1+u)
    (lambda^2*q2+lambda*q1)/q0
  }
  
  # cross-moment
  fmxm0=function(u){
    q0=Q(t,b0,u)
    q1=Q(t-1,b0,u)
    lambda=exp(b0+u)
    (q1/q0)^2*lambda^2
  }
  fmxm1=function(u){
    q0=Q(t,b0+b1,u)
    q1=Q(t-1,b0+b1,u)
    lambda=exp(b0+b1+u)
    (q1/q0)^2*lambda^2
  }
  
  # calculations
  set.seed(seed)
  S=1e6
  
  # control arm
  x=rnorm(S,0,sqrt(s0))
  mu0=mean(fmu0(x))
  mm0=mean(fmm0(x))
  tau0=mm0-(mu0)^2
  mxm0=mean(fmxm0(x))
  kappa0=sqrt(tau0)/mu0
  rho0=(mxm0-(mu0)^2)/tau0
  
  # intervention arm
  x=rnorm(S,0,sqrt(s1))
  mu1=mean(fmu1(x))
  mm1=mean(fmm1(x))
  tau1=mm1-(mu1)^2
  mxm1=mean(fmxm1(x))
  kappa1=sqrt(tau1)/mu1
  rho1=(mxm1-(mu1)^2)/tau1
  
  # effect size
  gamma0=log(mu0)
  gamma1=log(mu1/mu0)
  
  return(data.frame(gamma0=gamma0,gamma1=gamma1,mu0=mu0,mu1=mu1,
                    tau0=tau0,tau1=tau1,kappa0=kappa0,kappa1=kappa1,
                    rho0=rho0,rho1=rho1))
}

# example
calc(seed=0611,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=Inf)
calc(seed=0611,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=6)
calc(seed=0611,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=3)
calc(seed=0611,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=1)


calc(seed=0611,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=Inf)
calc(seed=0611,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=6)
calc(seed=0611,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=3)
calc(seed=0611,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=1)

# # check whether the results are stable
calc(seed=123,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=Inf)
calc(seed=123,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=6)
calc(seed=123,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=3)
calc(seed=123,b0=log(2.7),b1=log(0.7),s0=0.05,s1=0.1,t=1)
 


