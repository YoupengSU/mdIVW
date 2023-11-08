# Modified Debiased Inverse-Variance Weighted (mdIVW) Method for Mendelian Randomization

The modified debiased inverse-variance weighted (mdIVW) estimator is a Mendelian randomization method for estimating the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data. The mdIVW estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously.

## Setup
Use the following command in R to install the package (which is the latest version):
```
library(devtools)
install_github("YoupengSU/mdIVW",ref="main") 
```

## Usage
```
mr_mdIVW(
  gamma_hat,
  gamma_se,
  Gamma_hat,
  Gamma_se,
  rs_hat = NULL,
  rs_se = NULL,
  over.dispersion = TRUE,
  delta = 0,
  alpha = 0.05)
```
`gamma_hat`: A numeric vector of beta-coefficient values for genetic associations with the exposure variable.

`gamma_se`: The standard errors associated with the beta-coefficients `gamma_hat`.

`Gamma_hat`: A numeric vector of beta-coefficient values for genetic associations with the outcome variable.

`Gamma_se`: The standard errors associated with the beta-coefficients `Gamma_se`.	

`rs_hat`: A numeric vector of beta-coefficient values for genetic associations with the exposure variable, which will be used for the IV selection. rs_hat should be provided when delta is not zero.

`rs_se`: The standard errors associated with the beta-coefficients rs_hat.

`over.dispersion`: Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.

`delta`: The z-score threshold for IV selection. `delta` should be greater than or equal to zero. By default, `delta=0` (i.e., no IV selection will be conducted). 


`alpha`: The significance level used to calculate the confidence intervals. The default value is 0.05.

## Value
`Estimate`: The causal point estimate from the pIVW estimator.	

`StdError`: The standard error associated with `Estimate`.	

`CILower`: The lower bound of the confidence interval for `Estimate`.	

`CIUpper`: The upper bound of the confidence interval for `Estimate`. 	

`Pvalue`: P-value associated with `Estimate`.	

`Condition`: The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.	


## Example 
```
library(mr.mdIVW)  # load the mr.mdIVW package
attach(card_sol) # attach the real analysis data

# analyze the example data with the pIVW method. 
mr_mdIVW(gamma_hat = beta.exposure,
         gamma_se = se.exposure,
         Gamma_hat = beta.outcome,
         Gamma_se = se.outcome,
         rs_hat = rs_hat,
         rs_se = rs_se, 
         over.dispersion = TRUE,
         delta=0,
         alpha = 0.05)
detach(card_sol)  

# results 
Modified debiased inverse-variance weighted method

Over dispersion: TRUE 
IV selection threshold (delta): 0 
Number of variants : 2413

------------------------------------------------------------------
 Method Estimate Std Error 95% CI        p-value Condition
  mdIVW    0.982     0.243 0.506, 1.457 5.17e-05     5.054
------------------------------------------------------------------
```

