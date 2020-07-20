#####################################################################################
# Poisson SEE analysis of cluster randomized trials with truncated counts

# The following variance estimators are output from the program:
# (1) MB: Model-based estimator
# (2) BC0: Uncorrected sandwich estimator that extends Liang and Zeger, 1986
# (3) BC1: Bias-corrected sandwich estimator that extends Kauermann and Carroll, 2001
# (4) BC2: Bias-corrected sandwich estimator that extends Mancl and DeRouen, 2001
# (5) BC3: Bias-corrected sandwich estimator that extends Fay and Graubard, 2001

# June 2020

# INPUT
# y: The count outcome variable
# X: Marginal mean covariates (design matrix including intercept)
# id: Cluster identifier
# n: Vector of cluster sample sizes
# trt: cluster-level treatment indicator
# maxiter: Maximum number of iterations for Fisher Scoring
# epsilon: Tolerance for convergence

# Note that all inputs are required.
# ID's should be integers from 1 to K. Data should be sorted by ID before calling.
#####################################################################################

SEEPoisson=function(y, X, id, n, trt, maxiter=50, epsilon=0.00001){
  require(MASS)
  
  #####################################################################################
  # MODULE: BEGINEND
  # Creates two vectors that have the start and end points for each cluster
  
  # INPUT
  # n: Vector of cluster sample sizes
  
  # OUTPUT
  # first: Vector with starting row for cluster i
  # last: Vector with ending row for cluster i
  #####################################################################################
  
  BEGINEND=function(n){
    last=cumsum(n)
    first=last-n+1
    return(cbind(first,last))
  }
  
  #####################################################################################
  # Module: IS_POS_DEF
  # A = symmetric matrix
  # returns 1 if A is positive definite
  # 0 otherwise
  #####################################################################################
  
  is_pos_def=function(A){
    return(min(eigen(A)$values)>1e-13)
  }
  
  #####################################################################################
  # MODULE: SCORE
  # Generates the score matrix for each cluster and approximate information
  # to be used to estimate parameters and generate standard errors
  
  # INPUT
  # beta: Vector of marginal mean parameters
  # alpha: correlation vector
  # phi: dispersion vector
  # y: The count outcome
  # X: Marginal mean covariates
  # trt: cluster-level treatment indicator
  # n: Vector of cluster sample sizes
  # p: Number of marginal mean parameters
  
  # OUTPUT
  # U: Score vector
  # UUtran: Sum of U_i*U_i` across all clusters
  # Ustar: Approximate information matrix
  #####################################################################################
  
  SCORE=function(beta, alpha, phi, y, X, trt, n, p){
    U=rep(0,p)
    UUtran=Ustar=matrix(0,p,p)
    locx=BEGINEND(n)
    
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      
      U_c=rep(0,p)
      Ustar_c=matrix(0,p,p)
      mu_c=exp(c(X_c%*%beta))
      
      C=X_c*mu_c
      A=y_c-mu_c
      
      if(trt[i]==0){
        INVR=diag(1/(1-alpha[1]),n[i])-matrix(alpha[1]/((1-alpha[1])*(1-alpha[1]+n[i]*alpha[1])),n[i],n[i])
        INVB=diag(1,n[i]) %*% INVR %*% diag(1,n[i]) / phi[1]
      } else{
        INVR=diag(1/(1-alpha[2]),n[i])-matrix(alpha[2]/((1-alpha[2])*(1-alpha[2]+n[i]*alpha[2])),n[i],n[i])
        INVB=diag(1,n[i]) %*% INVR %*% diag(1,n[i]) / phi[2]
      }
      
      U_c=t(C)%*%INVB%*%A
      UUtran_c=tcrossprod(U_c)
      Ustar_c=t(C)%*%INVB%*%C
      U=U+U_c
      UUtran=UUtran+UUtran_c
      Ustar=Ustar+Ustar_c
    }
    return(list(U=U,UUtran=UUtran,Ustar=Ustar))
  }
  
  #####################################################################################
  # MODULE: INITBETA
  # Generates initial values for beta. 
  # independence poisson GEE
  
  # INPUT
  # y: The 0/1 binary outcomer4
  # X: Marginal mean covariates
  
  # OUTPUT
  # beta: Vector of marginal mean parameters
  #####################################################################################
  
  INITBETA=function(y,X){
    beta=as.numeric(coef(glm(y~-1+X,family=poisson(link=log))))
    return(beta)
  }
  
  # compute PHI
  getphi=function(y,X,trt,beta,n){
    n0=n[trt==0]
    n1=n[trt==1]
    locx=BEGINEND(n)
    phisum=phi=rep(0,2)
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      mu_c=exp(c(X_c%*%beta))
      rsq=(y_c-mu_c)^2
      if(trt[i]==0){
        phisum[1]=phisum[1]+sum(rsq)
      } else{
        phisum[2]=phisum[2]+sum(rsq)
      }
    }
    phi[1]=phisum[1]/sum(n0)
    phi[2]=phisum[2]/sum(n1)
    return(phi)
  }
  
  # compute ALPHA
  getalpha=function(y,X,trt,beta,phi,n){
    n0=n[trt==0]
    n1=n[trt==1]
    locx=BEGINEND(n)
    alpsum=alpha=rep(0,2)
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      mu_c=exp(c(X_c%*%beta))
      if(trt[i]==0){
        r=(y_c-mu_c)/sqrt(phi[1])
        rrtran=tcrossprod(r)
        alpsum[1]=alpsum[1]+sum(rrtran[upper.tri(rrtran)])
      } else{
        r=(y_c-mu_c)/sqrt(phi[2])
        rrtran=tcrossprod(r)
        alpsum[2]=alpsum[2]+sum(rrtran[upper.tri(rrtran)])
      }
    }
    alpha[1]=alpsum[1]/sum(n0*(n0-1)/2)
    alpha[2]=alpsum[2]/sum(n1*(n1-1)/2)
    return(alpha)
  }
  
  #####################################################################################
  # MODULE: INVBIG
  # compute (A - mm`)^{-1}c without performing the inverse directly
  # INPUT
  # ainvc: inverse of matrix A times vector c
  # ainvm: inverse of matrix A times matrix (with low no. of columns) M
  # M: matrix of eigen column vectors m1,m2, ..
  # c: vector 
  # start: of do loop
  # end:   of do loop, rank of X
  
  # OUTPUT
  # ainvc: inverse of matrix A times vector c
  #####################################################################################
  
  INVBIG=function(ainvc,ainvm,m,c,start,end){
    for(i in start:end){
      b=ainvm[,i]
      bt=t(b)
      btm=bt%*%m
      btmi=btm[,i]
      gam=1-btmi
      bg=b/gam
      ainvc=ainvc+bg%*%(bt%*%c)
      if(i<end){
        ainvm=ainvm+bg%*%btm
      }
    }
    return(ainvc)
  }
  
  #####################################################################################
  # MODULE: MAKEVAR
  # Creates covariance matrix of beta
  
  # INPUT
  # beta: Vector of marginal mean parameters
  # alpha: common correlation
  # phi: dispersion vector
  # y: The count outcome
  # X: Marginal mean covariates
  # trt: cluster-level treatment indicator
  # n: Vector of cluster sample sizes
  # p: Number of marginal mean parameters
  
  # OUTPUT
  # robust: Robust covariance matrix for beta and alpha
  # naive: Naive (Model-Based) covariance matrix for beta
  # varMD: bias-corrected variance by Mancl and Derouen (2001)
  # varKC: bias-corrected variance by Kauermann and Carroll (2001)
  # varFG: bias-corrected variance by Fay and Graubard (2001)
  # ROBFLAG: whether covariance is not positive definite 
  #####################################################################################
  
  MAKEVAR=function(beta, alpha, phi, y, X, trt, n, p){
    
    SCORE_RES=SCORE(beta, alpha, phi, y, X, trt, n, p)
    U=SCORE_RES$U
    UUtran=SCORE_RES$UUtran
    Ustar=SCORE_RES$Ustar
    
    naive=ginv(Ustar)
    
    # new commands to compute INV(I - H1)
    eigenRES1=eigen(naive)
    evals1=eigenRES1$values
    evecs1=eigenRES1$vectors
    sqrevals1=sqrt(evals1)
    sqe1=evecs1%*%diag(sqrevals1)
    
    # Bias-corrected variance
    Ustar_c_array=UUtran_c_array=array(0,c(p,p,length(n)))
    UUtran=UUbc=UUbc2=UUbc3=Ustar=matrix(0,p,p)
    
    locx=BEGINEND(n)
    
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      mu_c=exp(c(X_c%*%beta))
      
      U_i=U_c=rep(0,p)
      Ustar_c=matrix(0,p,p)
      
      # commands for beta
      C=X_c*mu_c
      A=y_c-mu_c
      if(trt[i]==0){
        INVR=diag(1/(1-alpha[1]),n[i])-matrix(alpha[1]/((1-alpha[1])*(1-alpha[1]+n[i]*alpha[1])),n[i],n[i])
        INVB=diag(1,n[i]) %*% INVR %*% diag(1,n[i]) / phi[1]
      } else{
        INVR=diag(1/(1-alpha[2]),n[i])-matrix(alpha[2]/((1-alpha[2])*(1-alpha[2]+n[i]*alpha[2])),n[i],n[i])
        INVB=diag(1,n[i]) %*% INVR %*% diag(1,n[i]) / phi[2]
      }
      U_i=t(C)%*%INVB%*%A
      
      # commands for generalized inverse - beta
      ai1=INVB
      mm1=C%*%sqe1
      ai1A=ai1%*%A
      ai1m1=ai1%*%mm1
      ai1A=INVBIG(ai1A,ai1m1,mm1,A,1,p)
      U_c=t(C)%*%ai1A
      
      Ustar_c=t(C)%*%INVB%*%C
      Ustar=Ustar+Ustar_c
      UUtran_c=tcrossprod(U_i)
      UUtran=UUtran+UUtran_c
      UUbc_c=tcrossprod(U_c)
      UUbc=UUbc+UUbc_c
      UUbc_ic=tcrossprod(U_c,U_i)
      UUbc2=UUbc2+UUbc_ic
      
      Ustar_c_array[,,i]=Ustar_c
      UUtran_c_array[,,i]=UUtran_c
    }
    
    # calculating adjustment factor for BC3
    for(i in 1:length(n)){      
      Hi=diag(1/sqrt(1-pmin(0.75,c(diag(Ustar_c_array[,,i]%*%naive)))))
      UUbc3=UUbc3+Hi%*%UUtran_c_array[,,i]%*%Hi
    }
    
    # BC0 or usual Sandwich estimator of Prentice (1988);     
    robust=naive%*%UUtran%*%naive
    
    # BC1 or Variance estimator due to Kauermann and Carroll (2001);
    varKC=naive%*%(UUbc2+t(UUbc2))%*%t(naive)/2
    
    # BC2 or Variance estimator due to Mancl and DeRouen (2001);
    varMD=naive%*%UUbc%*%t(naive)
    
    # BC3 or Variance estimator due to Fay and Graubard (2001);
    varFG=naive%*%UUbc3%*%t(naive)
    
    ROBFLAG=0
    if(min(diag(robust))<=0){ROBFLAG = 1}
    if(min(diag(varMD))<=0){ROBFLAG = 1}
    if(min(diag(varKC))<=0){ROBFLAG = 1}
    if(min(diag(varFG))<=0){ROBFLAG = 1}
    
    return(list(robust=robust,naive=naive,varMD=varMD,varKC=varKC,varFG=varFG,
                ROBFLAG=ROBFLAG))
  }
  
  #####################################################################################
  # MODULE: Performs modified Fisher scoring algorithm
  
  # INPUT
  # y: The count outcome
  # X: Marginal mean covariates
  # n: Vector of cluster sample sizes
  # maxiter: Max number of iterations
  # epsilon: Tolerence for convergence
  
  # OUTPUT
  # beta: A p x 1 vector of marginal mean parameter estimates  
  # alpha: A q x 1 vector of marginal correlation parameter estimates  
  # robust: Robust covariance matrix for beta and alpha
  # naive: Naive (Model-Based) covariance matrix for beta
  # varMD: Bias-corrected variance due to Mancl and DeRouen (2001)
  # varKC: Bias-corrected variance due to Kauermann and Carroll (2001)
  # varFG: Bias-corrected variance due to Fay and Graubard (2001)
  # niter: Number of iterations required for convergence
  # converge: Did the algorithm converge (0 = no, 1 = yes)
  #####################################################################################
  p=ncol(X)
  delta=rep(2*epsilon,p)
  deltaalp=deltaphi=rep(2*epsilon,2)
  converge=0
  beta=INITBETA(y,X)
  phi=getphi(y,X,trt,beta,n)
  alpha=getalpha(y,X,trt,beta,phi,n)
  
  niter=1
  while((niter<=maxiter) & (max(abs(c(delta,deltaalp,deltaphi)))>epsilon) ){
    SCORE_RES=SCORE(beta, alpha, phi, y, X, trt, n, p)
    U=SCORE_RES$U
    UUtran=SCORE_RES$UUtran
    Ustar=SCORE_RES$Ustar
    
    psdustar=is_pos_def(Ustar)
    mineig=min(eigen(Ustar)$values)
    if(psdustar==TRUE){
      delta=solve(Ustar,U)
      beta=beta+delta
    } else{
      SINGFLAG=1
    }
    phiold=phi
    phi=getphi(y,X,trt,beta,n)
    deltaphi=phi-phiold
    alphaold=alpha
    alpha=getalpha(y,X,trt,beta,phi,n)
    deltaalp=alpha-alphaold
    converge=(max(abs(c(delta,deltaalp,deltaphi)))<=epsilon)
    niter=niter+1
  }
  
  if(converge==1){
    MAKEVAR_RES=MAKEVAR(beta, alpha, phi, y, X, trt, n, p)
    robust=MAKEVAR_RES$robust
    naive=MAKEVAR_RES$naive
    varMD=MAKEVAR_RES$varMD
    varKC=MAKEVAR_RES$varKC
    varFG=MAKEVAR_RES$varFG
    ROBFLAG=MAKEVAR_RES$ROBFLAG
    return(list(beta=beta,phi=phi,alpha=alpha,naive=naive,robust=robust,varMD=varMD,
                varKC=varKC,varFG=varFG,niter=niter,converge=converge,
                ROBFLAG=ROBFLAG))
  } else{
    stop("Fisher scoring algorithm did not converge")
  }
}



