
##################################################
# simulation main code: empirical type I error or power
# 
# nsim: number of simulations
# mbar: mean cluster size
# b0: beta0 in Poisson mixed model (PMM)
# b1: beta1 in PMM
# s0: variance component in control arm
# s1: variance component in intervention arm
# t: truncation point
# cv: coefficient of variation of the cluster size distribution
# cor: correlation structure
#      "IND" = independence
#      "EXH" = arm-specific exchangeable 
###################################################

pw_simulation<-function(nsim, mbar, b0, b1 , s0, s1, t, cv, cor){
  library(geepack)
  
  # generate marginal parameters
  param<-calc(seed=0611,b0=b0,b1=b1,s0=s0,s1=s1,t=t)
  
  # estimate sample size
  ss<-sample_size(mu0=param[3], mu1=param[4], tau0=param[5], tau1=param[6], kappa0=param[7], kappa1=param[8],rho0=param[9], rho1=param[10], mbar=mbar, cv=cv, cor=cor)
  
  Err<-rep(0,nsim)
  output_alpha<-output_beta<-output_var_NV<-output_var_RB<-output_var_MD<-output_var_KC<-output_var_FG<-rep(NA,nsim)
  
  for(s in 1:nsim){
    if (s%%50==0) print(s)
    
    set.seed(100+s)
  
    df<-simdata(n=ss,mbar=mbar,cv=cv,b0=b0,b1=b1,s0=s0,s1=s1,t=t)
    
    cl_size<-as.numeric(table(df$datatest$cluster))
    error=0
    
    if (cor=="IND"){
      out=try(geeglm(df$datatest$y~-1+cbind(1,df$datatest$x),id=df$datatest$cluster,family=poisson(link="log"),corstr="independence"),silent=T)
    
    if(is(out, 'try-error')){
      Err[s]=1
      next
    } 
    
    res=summary(out)
    beta=as.numeric(coef(out))
    alpha=0
    phi=as.numeric(res$dispersion[1])
    fit<-LogPoisBCV(y=df$datatest$y,X=cbind(1,df$datatest$x),beta=beta,phi=phi,id=df$datatest$cluster)
    
    #output_alpha[s]=fit$alpha
    output_beta[s]=beta[2]
    output_var_NV[s]<-fit$naive[2,2]
    output_var_RB[s]<-fit$robust[2,2]
    output_var_MD[s]<-fit$varMD[2,2]
    output_var_KC[s]<-fit$varKC[2,2]
    output_var_FG[s]<-fit$varFG[2,2]    
    
    }
    
    if (cor=="EXH"){
      fit= try(SEEPoisson(y=df$datatest$y, X=cbind(1,df$datatest$x), id=df$datatest$cluster, n=cl_size, trt=df$trt, maxiter=50, epsilon=0.00001),silent=T)
    
      if(is(fit, 'try-error')){
        Err[s]=1
        next
      } 
      
    output_beta[s]=fit$beta[2]
    output_var_NV[s]<-fit$naive[2,2]
    output_var_RB[s]<-fit$robust[2,2]
    output_var_MD[s]<-fit$varMD[2,2]
    output_var_KC[s]<-fit$varKC[2,2]
    output_var_FG[s]<-fit$varFG[2,2]   
      
    }  
  }
  
  #alpha_mean<-mean(output_alpha,na.rm=T)
  
  test_stat_NV=output_beta/sqrt(output_var_NV)
  test_stat_RB=output_beta/sqrt(output_var_RB)
  test_stat_MD=output_beta/sqrt(output_var_MD)
  test_stat_KC=output_beta/sqrt(output_var_KC)
  test_stat_FG=output_beta/sqrt(output_var_FG)
  test_stat_MDKC=output_beta/(0.5*sqrt(output_var_KC)+0.5*sqrt(output_var_MD))
  
  # compute p-values
  PT_NV=2*(1-pt(abs(test_stat_NV),df=ss-2))
  PT_RB=2*(1-pt(abs(test_stat_RB),df=ss-2))
  PT_MD=2*(1-pt(abs(test_stat_MD),df=ss-2))
  PT_KC=2*(1-pt(abs(test_stat_KC),df=ss-2))
  PT_FG=2*(1-pt(abs(test_stat_FG),df=ss-2))
  PT_MDKC=2*(1-pt(abs(test_stat_MDKC),df=ss-2))
  
  pval=data.frame(PT_NV=PT_NV, PT_RB=PT_RB, PT_MD=PT_MD, PT_KC=PT_KC, PT_MDKC=PT_MDKC)
  Power<-cbind(c("PT_NV:","PT_RB:","PT_MD:","PT_KC:","PT_FG:","PT_MDKC:")
                 ,c(mean(PT_NV<0.05,na.rm=T),mean(PT_RB<0.05,na.rm=T),mean(PT_MD<0.05,na.rm=T),mean(PT_KC<0.05,na.rm=T),mean(PT_FG<0.05,na.rm=T),mean(PT_MDKC<0.05,na.rm=T)))
  
  colnames(Power)<-c("type","power")
  Power<-rbind(Power, c("err_rate",sum(Err)/length(Err)))
  
  print(c(nsim, mbar, s0, s1, t, cv, cor))
  print(param)
  print(Power)
  return(list(param=param, Power=Power))

}


##################################################
# simulation main code: empirical type I error
# 
# nsim: number of simulations
# mbar: mean cluster size
# b0: beta0 in Poisson mixed model (PMM)
# b1: beta1 in PMM
# s0: variance component in control arm
# s1: variance component in intervention arm
# t: truncation point
# cv: coefficient of variation of the cluster size distribution
# cor: correlation structure
#      "IND" = independence
#      "EXH" = arm-specific exchangeable 
###################################################

t1_simulation<-function(nsim, mbar, b0, b1, s0, s1, t, cv, cor){
  library(geepack)
  
  # generate marginal parameters
  param<-calc(seed=0611,b0=b0,b1=b1,s0=s0,s1=s1,t=t)
  
  # estimate sample size
  ss<-sample_size(mu0=param[3], mu1=param[4], tau0=param[5], tau1=param[6], kappa0=param[7], kappa1=param[8],rho0=param[9], rho1=param[10], mbar=mbar, cv=cv, cor=cor)
  
  # search for conditional RR parameter that satisfies marginal RR=0 (null hypothesis)
  opt_b1<-getb1(seed = 0611, b0=b0,s0=s0,s1=s1,t=t,step = 0.005)
  b1<-as.numeric(opt_b1$b1try)
  #param<-as.numeric(opt_b1$res)
  
  Err<-rep(0,nsim)
  output_alpha<-output_beta<-output_var_NV<-output_var_RB<-output_var_MD<-output_var_KC<-output_var_FG<-rep(NA,nsim)
  
  for(s in 1:nsim){
    if (s%%50==0) print(s)
    
    set.seed(100+s)
    
    df<-simdata(n=ss,mbar=mbar,cv=cv,b0=b0,b1=b1,s0=s0,s1=s1,t=t)
    
    cl_size<-as.numeric(table(df$datatest$cluster))
    error=0
    
    if (cor=="IND"){
      out=try(geeglm(df$datatest$y~-1+cbind(1,df$datatest$x),id=df$datatest$cluster,family=poisson(link="log"),corstr="independence"),silent=T)
      
      if(is(out, 'try-error')){
        Err[s]=1
        next
      } 
      
      res=summary(out)
      beta=as.numeric(coef(out))
      alpha=0
      phi=as.numeric(res$dispersion[1])
      fit<-LogPoisBCV(y=df$datatest$y,X=cbind(1,df$datatest$x),beta=beta,phi=phi,id=df$datatest$cluster)
      
      #output_alpha[s]=fit$alpha
      output_beta[s]=beta[2]
      output_var_NV[s]<-fit$naive[2,2]
      output_var_RB[s]<-fit$robust[2,2]
      output_var_MD[s]<-fit$varMD[2,2]
      output_var_KC[s]<-fit$varKC[2,2]
      output_var_FG[s]<-fit$varFG[2,2]    
      
    }
    
    if (cor=="EXH"){
      fit= try(SEEPoisson(y=df$datatest$y, X=cbind(1,df$datatest$x), id=df$datatest$cluster, n=cl_size, trt=df$trt, maxiter=50, epsilon=0.00001),silent=T)
      
      if(is(fit, 'try-error')){
        Err[s]=1
        next
      } 
      
      output_beta[s]=fit$beta[2]
      output_var_NV[s]<-fit$naive[2,2]
      output_var_RB[s]<-fit$robust[2,2]
      output_var_MD[s]<-fit$varMD[2,2]
      output_var_KC[s]<-fit$varKC[2,2]
      output_var_FG[s]<-fit$varFG[2,2]   
      
    }  
  }
  
  #alpha_mean<-mean(output_alpha,na.rm=T)
  
  test_stat_NV=output_beta/sqrt(output_var_NV)
  test_stat_RB=output_beta/sqrt(output_var_RB)
  test_stat_MD=output_beta/sqrt(output_var_MD)
  test_stat_KC=output_beta/sqrt(output_var_KC)
  test_stat_FG=output_beta/sqrt(output_var_FG)
  test_stat_MDKC=output_beta/(0.5*sqrt(output_var_KC)+0.5*sqrt(output_var_MD))
  
  # compute p-values
  PT_NV=2*(1-pt(abs(test_stat_NV),df=ss-2))
  PT_RB=2*(1-pt(abs(test_stat_RB),df=ss-2))
  PT_MD=2*(1-pt(abs(test_stat_MD),df=ss-2))
  PT_KC=2*(1-pt(abs(test_stat_KC),df=ss-2))
  PT_FG=2*(1-pt(abs(test_stat_FG),df=ss-2))
  PT_MDKC=2*(1-pt(abs(test_stat_MDKC),df=ss-2))
  
  
  pval=data.frame(PT_NV=PT_NV, PT_RB=PT_RB, PT_MD=PT_MD, PT_KC=PT_KC, PT_MDKC=PT_MDKC)
  t1<-cbind(c("PT_NV:","PT_RB:","PT_MD:","PT_KC:","PT_FG:","PT_MDKC:")
               ,c(mean(PT_NV<0.05,na.rm=T),mean(PT_RB<0.05,na.rm=T),mean(PT_MD<0.05,na.rm=T),mean(PT_KC<0.05,na.rm=T),mean(PT_FG<0.05,na.rm=T),mean(PT_MDKC<0.05,na.rm=T)))
  
  colnames(t1)<-c("type","t1")
  t1<-rbind(t1, c("err_rate",sum(Err)/length(Err)))
  
  print(c(nsim, mbar, s0, s1, t, cv, cor))
  print(param)
  print(t1)
  return(list(param=param, t1=t1))
  
}







