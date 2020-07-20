# TruncatedCounts

The current folder includes R code for data generation, sample size and power calculation, and implementing the GEE and SEE approaches in the article "Sample size and power considerations for cluster randomized trials with count outcomes subject to right truncation" by Li and Tong (under review).

(Truncated) count outcomes are generated from a conditional mean model with arm-specific variance components to mimic a parallel cluster randomized trial. The treatment is assiged at the cluster level.

List of Files:
1) simdata.r = R file for simulating correlated count data subject to right truncation;
2) calc.r = R file for implementing the set of conversion formulas to obtain marginal quantities from the conditional model parameters. These marginal quantities are input parameters for sample size calculation. The details of the conversion formulas are presented in the Appendix of Li and Tong (under review).
3) sample_size_calculation.r = R file for implementing the sample size formula developed under the independence working correlation and the arm-specific exchangeable working correlation. The sample size is rounded to the nearest even integer above.
4) LogPoisBCV.r = R file for calculating the bias-corrected sandwich variances for the treatment effect under the Poisson GEE with the independence working correlation model.
5) SEEPoisson.r = R file for implementing the SEE approach and calculating the bias-corrected sandwich variances for the treatment effect. The working correlation model is the arm-specific exchangeable structure. 
