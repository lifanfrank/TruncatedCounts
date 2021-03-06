# TruncatedCounts

The current folder includes R code for reproducing Table 1-5 and Figure 1 in the article "Sample size and power considerations for cluster randomized trials with count outcomes subject to right truncation" by Li and Tong (Biometrical Journal, 2021).

For questions or comments about the code please contact Fan Li at <fan.f.li@yale.edu>.

I. List of Supporting Files: These supporting files are sourced in the main files that reproduce the numbers in the submitted manuscript

1) simdata.r = function to simulate correlated count data subject to right truncation;
2) calc.r = function to implement the set of conversion formulas to obtain marginal quantities from the conditional model parameters. These marginal quantities are input parameters for sample size calculation. The details of the conversion formulas are presented in the Appendix;
3) optimal_b1.r = function to numerically determine the conditional RR parameter that leads to zero marginal RR. This function is used to determine the simulation null parameters;
4) sample_size_calculation.r = function to implement the sample size formula developed under the independence working correlation and the arm-specific exchangeable working correlation. The sample size is rounded to the nearest even integer above;
5) sample_size_calculation_tp.r = same as 4), except that this function also produces the exact predicted power from the formula;
6) LogPoisBCV.r = function to calculate the bias-corrected sandwich variances for the treatment effect under the Poisson GEE with the independence working correlation model;
7) SEEPoisson.r = function to implement the SEE approach and calculate the bias-corrected sandwich variances for the treatment effect. The working correlation model is the arm-specific exchangeable model;
8) power_t1_calculation.r = function to execute the simulation study and compute the empirical type I error rate or power. 

II. List of Main Files: These main files are used to reproduce the results in Table 1-5 and Figure 1.

9) res1_Simulation_Table1_and_Table2.r = reproduce simulation results in Table 1 and Table 2;
10) res2_Application_Table3_to_Table5.r = reproduce numerical results in Table 3 to Table 5;
11) res3_Application_Figure1.r = reproduce numbers and the format of Figure 1.

III. Software 

Analyses were conducted with R, version 3.6.1 (https://www.r-project.org/)
The calculations used R packages geepack (version 1.2-1), truncdist (version 1.0-2), MASS (version 7.3-51.4), xtable (version 1.8-4) and plotrix (version 3.7-8)

IV. R commands for the installation of R packages 

install.packages(c("geepack", "truncdist", "MASS", "xtable", "plotrix")) 
