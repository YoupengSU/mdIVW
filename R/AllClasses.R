#' PIVW Class
#'
#' @description An object containing the estimate produced using the penalized inverse-variance weighted (pIVW) method as well as various statistics.
#'
#' @slot Over.dispersion Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.
#' @slot Delta The z-score threshold for IV selection. By default, \code{delta=0} (i.e., no IV selection will be conducted).
#' @slot Estimate The causal point estimate from the pIVW estimator.
#' @slot StdError The standard error associated with \code{Estimate}.
#' @slot CILower The lower bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then lower limits of all ranges will be reported.
#' @slot CIUpper The upper bound of the confidence interval for \code{Estimate}, which is derived from the bootstrapping Fieller method or normal distribution. For the bootstrapping Fieller's interval, if it contains multiple ranges, then upper limits of all ranges will be reported.
#' @slot Alpha The significance level used in constructing the confidence interval (default is 0.05).
#' @slot Pvalue P-value associated with the causal estimate from the pIVW estimator.
#' @slot Tau2 The variance of the balanced horizontal pleiotropy. \code{Tau2} is calculated by using all IVs in the data before conducting the IV selection.
#' @slot SNPs The number of SNPs after IV selection.
#' @slot Condition The estimated effective sample size. It is recommended to be greater than 5 for the pIVW estimator to achieve reliable asymptotic properties.
#'
#' @keywords internal


setClass("mdIVW",
         representation(Over.dispersion = "logical",


                        Delta = "numeric",

                        Estimate = "numeric",
                        StdError = "numeric",
                        CILower = "numeric",
                        CIUpper = "numeric",
                        Alpha = "numeric",

                        Pvalue = "numeric",
                        Tau2 = "numeric",
                        SNPs = "numeric",
                        Condition = "numeric")
)















