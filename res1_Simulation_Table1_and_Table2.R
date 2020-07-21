#################################################
# This R program sources all the components in 
# the folder and reproduces Table 1 and Table 2 
# of the manuscript
#################################################

setwd(".../ReproducibleResearch")  # need to change the directory here

source("calc.R")
source("simdata.R")
source("LogPoisBCV.R")
source("optimal_b1.R")
source("power_t1_calculation.R")
source("sample_size_calculation.R")
source("SEEPoisson.R")

require(geepack)
require(truncdist)
require(MASS)
require(xtable)

########################################################################
# The following code reproduces Table 1 on empirical type I error rate

# PT_NV: t-test with model-based variance (MB)
# PT_RB: t-test with Liang and Zeger variance (LZ)
# PT_MD: t-test with Mancl and DeRouen variance (MD)
# PT_KC: t-test with Kauermann and Carroll variance (KC)
# PT_MDKC: t-test with the average MD/KC standard error (AVG)
# err_rate: number of non-convergences
########################################################################

#########################
# Left column of Table 1
#########################

# Block 1: Independence working correlation (s0=0.05, s1=0.05)
T1_M25_B125055_S0505_T9_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=Inf, cv=0.0, cor="IND")
T1_M25_B125055_S0505_T4_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=4  , cv=0.0, cor="IND")
T1_M25_B125055_S0505_T2_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=2  , cv=0.0, cor="IND")
T1_M25_B125055_S0505_T1_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=1  , cv=0.0, cor="IND")

# Block 2: Independence working correlation (s0=0.05, s1=0.10)
T1_M25_B125055_S0510_T9_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=Inf, cv=0.0, cor="IND")
T1_M25_B125055_S0510_T4_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=4  , cv=0.0, cor="IND")
T1_M25_B125055_S0510_T2_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=2  , cv=0.0, cor="IND")
T1_M25_B125055_S0510_T1_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=1  , cv=0.0, cor="IND")

# Block 3: Independence working correlation (s0=0.05, s1=0.20)
T1_M25_B125055_S0520_T9_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=Inf, cv=0.0, cor="IND")
T1_M25_B125055_S0520_T4_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=4  , cv=0.0, cor="IND")
T1_M25_B125055_S0520_T2_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=2  , cv=0.0, cor="IND")
T1_M25_B125055_S0520_T1_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=1  , cv=0.0, cor="IND")

# Block 4: Independence working correlation (s0=0.10, s1=0.10)
T1_M25_B125055_S1010_T9_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=Inf, cv=0.0, cor="IND")
T1_M25_B125055_S1010_T4_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=4  , cv=0.0, cor="IND")
T1_M25_B125055_S1010_T2_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=2  , cv=0.0, cor="IND")
T1_M25_B125055_S1010_T1_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=1  , cv=0.0, cor="IND")

# Block 5: Independence working correlation (s0=0.10, s1=0.20)
T1_M25_B125055_S1020_T9_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=Inf, cv=0.0, cor="IND")
T1_M25_B125055_S1020_T4_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=4  , cv=0.0, cor="IND")
T1_M25_B125055_S1020_T2_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=2  , cv=0.0, cor="IND")
T1_M25_B125055_S1020_T1_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=1  , cv=0.0, cor="IND")

# Block 6: Independence working correlation (s0=0.20, s1=0.20)
T1_M25_B125055_S2020_T9_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=Inf, cv=0.0, cor="IND")
T1_M25_B125055_S2020_T4_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=4  , cv=0.0, cor="IND")
T1_M25_B125055_S2020_T2_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=2  , cv=0.0, cor="IND")
T1_M25_B125055_S2020_T1_CV00_IND<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=1  , cv=0.0, cor="IND")


#########################
# Right column of Table 1
#########################

# Block 1: arm-specific exchangeable correlation (s0=0.05, s1=0.05)
T1_M25_B125055_S0505_T9_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=Inf, cv=0.0, cor="EXH")
T1_M25_B125055_S0505_T4_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=4  , cv=0.0, cor="EXH")
T1_M25_B125055_S0505_T2_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=2  , cv=0.0, cor="EXH")
T1_M25_B125055_S0505_T1_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=1  , cv=0.0, cor="EXH")

# Block 2: arm-specific exchangeable correlation (s0=0.05, s1=0.10)
T1_M25_B125055_S0510_T9_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=Inf, cv=0.0, cor="EXH")
T1_M25_B125055_S0510_T4_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=4  , cv=0.0, cor="EXH")
T1_M25_B125055_S0510_T2_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=2  , cv=0.0, cor="EXH")
T1_M25_B125055_S0510_T1_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=1  , cv=0.0, cor="EXH")

# Block 3: arm-specific exchangeable correlation (s0=0.05, s1=0.20)
T1_M25_B125055_S0520_T9_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=Inf, cv=0.0, cor="EXH")
T1_M25_B125055_S0520_T4_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=4  , cv=0.0, cor="EXH")
T1_M25_B125055_S0520_T2_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=2  , cv=0.0, cor="EXH")
T1_M25_B125055_S0520_T1_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=1  , cv=0.0, cor="EXH")

# Block 4: arm-specific exchangeable correlation (s0=0.10, s1=0.10)
T1_M25_B125055_S1010_T9_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=Inf, cv=0.0, cor="EXH")
T1_M25_B125055_S1010_T4_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=4  , cv=0.0, cor="EXH")
T1_M25_B125055_S1010_T2_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=2  , cv=0.0, cor="EXH")
T1_M25_B125055_S1010_T1_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=1  , cv=0.0, cor="EXH")

# Block 5: arm-specific exchangeable correlation (s0=0.10, s1=0.20)
T1_M25_B125055_S1020_T9_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=Inf, cv=0.0, cor="EXH")
T1_M25_B125055_S1020_T4_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=4  , cv=0.0, cor="EXH")
T1_M25_B125055_S1020_T2_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=2  , cv=0.0, cor="EXH")
T1_M25_B125055_S1020_T1_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=1  , cv=0.0, cor="EXH")

# Block 6: arm-specific exchangeable correlation (s0=0.20, s1=0.20)
T1_M25_B125055_S2020_T9_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=Inf, cv=0.0, cor="EXH")
T1_M25_B125055_S2020_T4_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=4  , cv=0.0, cor="EXH")
T1_M25_B125055_S2020_T2_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=2  , cv=0.0, cor="EXH")
T1_M25_B125055_S2020_T1_CV00_EXH<-t1_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=1  , cv=0.0, cor="EXH")





########################################################################
# The following code reproduces Table 2 on empirical power

# For this part, we need to calculate the predicted power (by formula) under each scenario
source("sample_size_calculation_tp.R")

# collection of right truncation points
t1<-c(Inf,4,2,1)

# Generate a table with estimated sample sizes and predicted power (by formula)
ss_table<-function(TT=t1,CV=0){
  
  samplesize<-matrix(NA,6,8)
  truepower<-matrix(NA,6,8)
  
  # set of variance components
  s<-matrix(c(0.05,0.05,0.05,0.10,0.05,0.20,0.10,0.10,0.10,0.20,0.20,0.20),2,6)
  corr<-c("IND","EXH")
  
  for (i in 1:6){
    temp1<-c()
    temp2<-c()
    for (corr in c("IND","EXH")){
      for (t in TT){
        ccc<-calc(seed=0611,b0=log(1.25), b1=log(0.55), s0=s[1,i], s1=s[2,i],t=t)
        ss<-sample_size_tp(ccc[3], ccc[4], ccc[5], ccc[6], ccc[7], ccc[8],ccc[9], ccc[10], mbar=25, cv=CV, cor=corr)  
        temp1<-c(temp1,ss[1])
        temp2<-c(temp2,ss[2])
      }
    }
    samplesize[i,]<-temp1
    truepower[i,]<-temp2
  }  
  return(list(samplesize,truepower))
}

# Table of sample size and predicted power
ssA00<-ss_table(T=t1,CV=0)

#######################################################################
# Compute the difference between empirical and predicted power in Table 2

# PT_NV: t-test with model-based variance (MB)
# PT_RB: t-test with Liang and Zeger variance (LZ)
# PT_MD: t-test with Mancl and DeRouen variance (MD)
# PT_KC: t-test with Kauermann and Carroll variance (KC)
# PT_MDKC: t-test with the average MD/KC standard error (AVG)
# err_rate: number of non-convergences
########################################################################


############################################
# Left column of Table 2 (empirical power)
############################################

# Block 1: Independence working correlation (s0=0.05, s1=0.05)
PW_M25_B125055_S0505_T9_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=Inf, cv=0.0, cor="IND")
PW_M25_B125055_S0505_T4_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=4  , cv=0.0, cor="IND")
PW_M25_B125055_S0505_T2_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=2  , cv=0.0, cor="IND")
PW_M25_B125055_S0505_T1_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=1  , cv=0.0, cor="IND")

# Block 2: Independence working correlation (s0=0.05, s1=0.10)
PW_M25_B125055_S0510_T9_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=Inf, cv=0.0, cor="IND")
PW_M25_B125055_S0510_T4_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=4  , cv=0.0, cor="IND")
PW_M25_B125055_S0510_T2_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=2  , cv=0.0, cor="IND")
PW_M25_B125055_S0510_T1_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=1  , cv=0.0, cor="IND")

# Block 3: Independence working correlation (s0=0.05, s1=0.20)
PW_M25_B125055_S0520_T9_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=Inf, cv=0.0, cor="IND")
PW_M25_B125055_S0520_T4_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=4  , cv=0.0, cor="IND")
PW_M25_B125055_S0520_T2_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=2  , cv=0.0, cor="IND")
PW_M25_B125055_S0520_T1_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=1  , cv=0.0, cor="IND")

# Block 4: Independence working correlation (s0=0.10, s1=0.10)
PW_M25_B125055_S1010_T9_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=Inf, cv=0.0, cor="IND")
PW_M25_B125055_S1010_T4_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=4  , cv=0.0, cor="IND")
PW_M25_B125055_S1010_T2_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=2  , cv=0.0, cor="IND")
PW_M25_B125055_S1010_T1_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=1  , cv=0.0, cor="IND")

# Block 5: Independence working correlation (s0=0.10, s1=0.20)
PW_M25_B125055_S1020_T9_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=Inf, cv=0.0, cor="IND")
PW_M25_B125055_S1020_T4_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=4  , cv=0.0, cor="IND")
PW_M25_B125055_S1020_T2_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=2  , cv=0.0, cor="IND")
PW_M25_B125055_S1020_T1_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=1  , cv=0.0, cor="IND")

# Block 6: Independence working correlation (s0=0.20, s1=0.20)
PW_M25_B125055_S2020_T9_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=Inf, cv=0.0, cor="IND")
PW_M25_B125055_S2020_T4_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=4  , cv=0.0, cor="IND")
PW_M25_B125055_S2020_T2_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=2  , cv=0.0, cor="IND")
PW_M25_B125055_S2020_T1_CV00_IND<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=1  , cv=0.0, cor="IND")



############################################
# Right column of Table 2 (empirical power)
############################################

# Block 1: Independence working correlation (s0=0.05, s1=0.05)
PW_M25_B125055_S0505_T9_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=Inf, cv=0.0, cor="EXH")
PW_M25_B125055_S0505_T4_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=4  , cv=0.0, cor="EXH")
PW_M25_B125055_S0505_T2_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=2  , cv=0.0, cor="EXH")
PW_M25_B125055_S0505_T1_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.05, t=1  , cv=0.0, cor="EXH")

# Block 2: Independence working correlation (s0=0.05, s1=0.10)
PW_M25_B125055_S0510_T9_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=Inf, cv=0.0, cor="EXH")
PW_M25_B125055_S0510_T4_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=4  , cv=0.0, cor="EXH")
PW_M25_B125055_S0510_T2_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=2  , cv=0.0, cor="EXH")
PW_M25_B125055_S0510_T1_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.10, t=1  , cv=0.0, cor="EXH")

# Block 3: Independence working correlation (s0=0.05, s1=0.20)
PW_M25_B125055_S0520_T9_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=Inf, cv=0.0, cor="EXH")
PW_M25_B125055_S0520_T4_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=4  , cv=0.0, cor="EXH")
PW_M25_B125055_S0520_T2_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=2  , cv=0.0, cor="EXH")
PW_M25_B125055_S0520_T1_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.05, s1=0.20, t=1  , cv=0.0, cor="EXH")

# Block 4: Independence working correlation (s0=0.10, s1=0.10)
PW_M25_B125055_S1010_T9_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=Inf, cv=0.0, cor="EXH")
PW_M25_B125055_S1010_T4_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=4  , cv=0.0, cor="EXH")
PW_M25_B125055_S1010_T2_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=2  , cv=0.0, cor="EXH")
PW_M25_B125055_S1010_T1_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.10, t=1  , cv=0.0, cor="EXH")

# Block 5: Independence working correlation (s0=0.10, s1=0.20)
PW_M25_B125055_S1020_T9_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=Inf, cv=0.0, cor="EXH")
PW_M25_B125055_S1020_T4_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=4  , cv=0.0, cor="EXH")
PW_M25_B125055_S1020_T2_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=2  , cv=0.0, cor="EXH")
PW_M25_B125055_S1020_T1_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.10, s1=0.20, t=1  , cv=0.0, cor="EXH")

# Block 6: Independence working correlation (s0=0.20, s1=0.20)
PW_M25_B125055_S2020_T9_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=Inf, cv=0.0, cor="EXH")
PW_M25_B125055_S2020_T4_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=4  , cv=0.0, cor="EXH")
PW_M25_B125055_S2020_T2_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=2  , cv=0.0, cor="EXH")
PW_M25_B125055_S2020_T1_CV00_EXH<-pw_simulation(1000, mbar=25, b0=log(1.25), b1=log(0.55), s0=0.20, s1=0.20, t=1  , cv=0.0, cor="EXH")

###################################################################
# Summarizing difference between empirical and predicted power
###################################################################

# write a function to summarize the results in Table 1 and Table 2
output<-function(PT="T1",B="125055",S="0505",CV="00"){
  require(xtable)  
  
  # compute predicted power * 100
  tp<-ssA00[[2]]*100
  
  if (S=="0505") tp0<-tp[1,]
  if (S=="0510") tp0<-tp[2,]
  if (S=="0520") tp0<-tp[3,]
  if (S=="1010") tp0<-tp[4,]
  if (S=="1020") tp0<-tp[5,]
  if (S=="2020") tp0<-tp[6,]
  
  df1<-get(paste(PT,"_M25_B",B,"_S",S,"_T9_CV",CV,"_IND",sep=""))
  df2<-get(paste(PT,"_M25_B",B,"_S",S,"_T4_CV",CV,"_IND",sep=""))
  df3<-get(paste(PT,"_M25_B",B,"_S",S,"_T2_CV",CV,"_IND",sep=""))
  df4<-get(paste(PT,"_M25_B",B,"_S",S,"_T1_CV",CV,"_IND",sep=""))
  df5<-get(paste(PT,"_M25_B",B,"_S",S,"_T9_CV",CV,"_EXH",sep=""))
  df6<-get(paste(PT,"_M25_B",B,"_S",S,"_T4_CV",CV,"_EXH",sep=""))
  df7<-get(paste(PT,"_M25_B",B,"_S",S,"_T2_CV",CV,"_EXH",sep=""))
  df8<-get(paste(PT,"_M25_B",B,"_S",S,"_T1_CV",CV,"_EXH",sep=""))
  
  names<-c("MB","LZ","MD","KC","FG","AVG")
  
  if (PT=="T1"){
    
    tab1<-cbind(round(as.numeric(df1$t1[1:6,2])*100,1),round(as.numeric(df2$t1[1:6,2])*100,1),round(as.numeric(df3$t1[1:6,2])*100,1),round(as.numeric(df4$t1[1:6,2])*100,1),
                round(as.numeric(df5$t1[1:6,2])*100,1),round(as.numeric(df6$t1[1:6,2])*100,1),round(as.numeric(df7$t1[1:6,2])*100,1),round(as.numeric(df8$t1[1:6,2])*100,1))
  }else{
    
    tab1<-cbind(round(as.numeric(df1$Power[1:6,2])*100,1),round(as.numeric(df2$Power[1:6,2])*100,1),round(as.numeric(df3$Power[1:6,2])*100,1),round(as.numeric(df4$Power[1:6,2])*100,1),
                round(as.numeric(df5$Power[1:6,2])*100,1),round(as.numeric(df6$Power[1:6,2])*100,1),round(as.numeric(df7$Power[1:6,2])*100,1),round(as.numeric(df8$Power[1:6,2])*100,1))
    
    #tab0=tab1
    for (i in 1:dim(tab1)[1]) tab1[i,]<-round(tab1[i,]-tp0,1)
    
  }
  
  #print(tab0)
  tab1<-data.frame(X0=NA,names,tab1,row.names=NULL)
  print(xtable(tab1,digits = 1), include.rownames = FALSE)
}  

# Type I error rates 
# (Latex code for Table 1)

output(PT="T1",B="125055",S="0505",CV="00")
output(PT="T1",B="125055",S="0510",CV="00")
output(PT="T1",B="125055",S="0520",CV="00")
output(PT="T1",B="125055",S="1010",CV="00")
output(PT="T1",B="125055",S="1020",CV="00")
output(PT="T1",B="125055",S="2020",CV="00")

# Difference between empirical and predicted power 
# (Latex code for Table 2)

output(PT="PW",B="125055", S="0505",CV="00")
output(PT="PW",B="125055", S="0510",CV="00")
output(PT="PW",B="125055", S="0520",CV="00")
output(PT="PW",B="125055", S="1010",CV="00")
output(PT="PW",B="125055", S="1020",CV="00")
output(PT="PW",B="125055", S="2020",CV="00")

