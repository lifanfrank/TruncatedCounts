#####################################################
# numerically determine the conditional RR parameter
# that leads to zero marginal RR
# 
# this function is used to determine the simulation
# null parameters
######################################################


getb1=function(seed=NULL,b0,s0,s1,t=Inf,step=0.005,epsilon=0.0005){
  b1new=b1old=log(1)
  delta=1
  resold=calc(seed,b0,b1=b1old,s0,s1,t)
  gamma1old=resold$gamma1
  
  while(delta>=0){
    b1old=b1new
    if(gamma1old>0){
      b1new=b1old-step
    } else{
      b1new=b1old+step
    }
    resnew=calc(seed,b0,b1=b1new,s0,s1,t)
    gamma1new=resnew$gamma1
    delta=gamma1old*gamma1new
    gamma1old=gamma1new
  }
  bound=sort(c(b1old,b1new))
  lb=bound[1]
  ub=bound[2]
  i=0
  b1try=b1new
  while(abs(gamma1new)>=epsilon){
    b1try=(lb+ub)/2
    restry=calc(seed,b0,b1=b1try,s0,s1,t)
    gamma1new=restry$gamma1
    if(gamma1new>0){
      ub=b1try
    } else{
      lb=b1try
    }
    i=i+1
    print(i)
  }
  res=calc(seed,b0,b1=b1try,s0,s1,t)
  return(list(b1try=b1try,res=res))
}