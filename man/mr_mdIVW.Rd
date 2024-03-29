% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr.mdIVW.R
\name{mr_mdIVW}
\alias{mr_mdIVW}
\title{modified debiasd Inverse-Variance Weighted (mdIVW) Method for Mendelian Randomization}
\usage{
mr_mdIVW(
  gamma_hat,
  gamma_se,
  Gamma_hat,
  Gamma_se,
  rs_hat = NULL,
  rs_se = NULL,
  over.dispersion = TRUE,
  delta = 0,
  alpha = 0.05
)
}
\arguments{
\item{gamma_hat}{A numeric vector of beta-coefficient values for genetic associations with the exposure variable.}

\item{gamma_se}{The standard errors associated with the beta-coefficients \code{gamma_hat}.}

\item{Gamma_hat}{A numeric vector of beta-coefficient values for genetic associations with the outcome variable.}

\item{Gamma_se}{The standard errors associated with the beta-coefficients \code{Gamma_hat}.}

\item{rs_hat}{A numeric vector of beta-coefficient values for genetic associations with the exposure variable, which will be used for the IV selection. \code{rs_hat} should be provided when \code{delta} is not zero. See 'Details'.}

\item{rs_se}{The standard errors associated with the beta-coefficients \code{rs_hat}.}

\item{over.dispersion}{Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.}

\item{delta}{The z-score threshold for IV selection. \code{delta} should be greater than or equal to zero. By default, \code{delta=0} (i.e., no IV selection will be conducted).  See 'Details'.}

\item{alpha}{The significance level used to calculate the confidence intervals. The default value is 0.05.}
}
\value{
The output from the function is a \code{mdIVW} object containing:

 \item{Over.dispersion}{\code{TRUE} if the method has considered balanced horizontal pleiotropy, \code{FALSE} otherwise.}
 \item{Delta}{The z-score threshold for IV selection.}
 \item{Estimate}{The causal point estimate from the mdIVW estimator.}
 \item{StdError}{The standard error associated with \code{Estimate}.}
 \item{CILower}{The lower bound of the confidence interval for \code{Estimate}, which is derived from the normal distribution.}
 \item{CIUpper}{The upper bound of the confidence interval for \code{Estimate}, which is derived from the normal distribution.}
 \item{Alpha}{The significance level used in constructing the confidence interval.}
 \item{Pvalue}{P-value associated with \code{Estimate}, which is derived from the normal distribution.}
 \item{Tau2}{The variance of the balanced horizontal pleiotropy. \code{Tau2} is calculated By using all IVs in the data before conducting the IV selection.}
 \item{SNPs}{The number of SNPs after IV selection.}
 \item{Condition}{The estimated effective sample size. It is recommended to be greater than 10 for the mdIVW estimator to achieve reliable asymptotic properties. See 'Details'.}
}
\description{
The modified debiasd Inverse-Variance Weighted (mdIVW) estimator is a Mendelian randomization method for estimating the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data. The mdIVW estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously.
}
\details{
The modified debiasd Inverse-Variance Weighted (mdIVW) estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously
in two-sample MR with summary statistics, i.e., an exposure sample (with IV-exposure effect \code{gamma_hat} and standard error \code{gamma_se}) and
an outcome sample (with IV-outcome effect \code{Gamma_hat} and standard error \code{Gamma_se}).

The mdIVW estimator also allows for IV selection in three-sample MR, where weak IVs are screened out using
an extra sample (with IV-exposure effect \code{rs_hat} and standard error \code{rs_se}) independent of the exposure sample and outcome sample.
Generally, the P-value for \code{rs_hat} can be computed By \code{sel.pval=2*pnorm(abs(rs_hat/rs_se), lower.tail = FALSE)},
Given \code{sel.pval} and a z-score threshold \code{delta}, the variants kept in the analysis will be those
with \code{sel.pval<2*pnorm(delta,lower.tail = FALSE)}.

The \code{mr_mdivw} function outputs a measure \code{Condition} that needs to be large for reliable asymptotic properties of the mdIVW estimator.
We also refer to \code{Condition} as effective sample size, which is a function of a measure of IV strength and the number of IVs.
When \code{delta} is zero (i.e., no IV selection), \code{Condition = (average F-statistic -1)*sqrt(# snps)}. When \code{delta} is not zero
(i.e., IV selection is conducted), \code{Condition = [(average F-statistic -1)*sqrt(# snps)]/phi_hat},
where the numerator is computed using the selected variants, and the denominator \code{phi_hat} involves the selection probabilities
of all variants. We suggest that \code{Condition} should be greater than 10 for the mdIVW estimator to achieve reliable asymptotic properties.
}
\examples{
mr_mdIVW(gamma_hat = Bx_exp, gamma_se = Bxse_exp, Gamma_hat = By_exp, Gamma_se = Byse_exp)

}
