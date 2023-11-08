#' modified debiasd Inverse-Variance Weighted (mdIVW) Method for Mendelian Randomization
#'

#' @description The modified debiasd Inverse-Variance Weighted (mdIVW) estimator is a Mendelian randomization method for estimating the causal effect of an exposure variable on an outcome of interest based on summary-level GWAS data. The mdIVW estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously.
#'
#' @param gamma_hat A numeric vector of beta-coefficient values for genetic associations with the exposure variable.
#' @param gamma_se The standard errors associated with the beta-coefficients \code{gamma_hat}.
#' @param Gamma_hat A numeric vector of beta-coefficient values for genetic associations with the outcome variable.
#' @param Gamma_se The standard errors associated with the beta-coefficients \code{Gamma_hat}.
#' @param rs_hat A numeric vector of beta-coefficient values for genetic associations with the exposure variable, which will be used for the IV selection. \code{rs_hat} should be provided when \code{delta} is not zero. See 'Details'.
#' @param rs_se  The standard errors associated with the beta-coefficients \code{rs_hat}.
#' @param over.dispersion Should the method consider overdispersion (balanced horizontal pleiotropy)? Default is TRUE.
#' @param delta The z-score threshold for IV selection. \code{delta} should be greater than or equal to zero. By default, \code{delta=0} (i.e., no IV selection will be conducted).  See 'Details'.
#' @param alpha The significance level used to calculate the confidence intervals. The default value is 0.05.
#'
#' @details The modified debiasd Inverse-Variance Weighted (mdIVW) estimator accounts for weak instruments and balanced horizontal pleiotropy simultaneously
#' in two-sample MR with summary statistics, i.e., an exposure sample (with IV-exposure effect \code{gamma_hat} and standard error \code{gamma_se}) and
#' an outcome sample (with IV-outcome effect \code{Gamma_hat} and standard error \code{Gamma_se}).
#'
#' The mdIVW estimator also allows for IV selection in three-sample MR, where weak IVs are screened out using
#' an extra sample (with IV-exposure effect \code{rs_hat} and standard error \code{rs_se}) independent of the exposure sample and outcome sample.
#' Generally, the P-value for \code{rs_hat} can be computed By \code{sel.pval=2*pnorm(abs(rs_hat/rs_se), lower.tail = FALSE)},
#' Given \code{sel.pval} and a z-score threshold \code{delta}, the variants kept in the analysis will be those
#' with \code{sel.pval<2*pnorm(delta,lower.tail = FALSE)}.
#'
#' The \code{mr_mdivw} function outputs a measure \code{Condition} that needs to be large for reliable asymptotic properties of the mdIVW estimator.
#' We also refer to \code{Condition} as effective sample size, which is a function of a measure of IV strength and the number of IVs.
#' When \code{delta} is zero (i.e., no IV selection), \code{Condition = (average F-statistic -1)*sqrt(# snps)}. When \code{delta} is not zero
#' (i.e., IV selection is conducted), \code{Condition = [(average F-statistic -1)*sqrt(# snps)]/phi_hat},
#' where the numerator is computed using the selected variants, and the denominator \code{phi_hat} involves the selection probabilities
#' of all variants. We suggest that \code{Condition} should be greater than 10 for the mdIVW estimator to achieve reliable asymptotic properties.
#'
#'
#' @return The output from the function is a \code{mdIVW} object containing:
#'
#'  \item{Over.dispersion}{\code{TRUE} if the method has considered balanced horizontal pleiotropy, \code{FALSE} otherwise.}
#'  \item{Delta}{The z-score threshold for IV selection.}
#'  \item{Estimate}{The causal point estimate from the mdIVW estimator.}
#'  \item{StdError}{The standard error associated with \code{Estimate}.}
#'  \item{CILower}{The lower bound of the confidence interval for \code{Estimate}, which is derived from the normal distribution.}
#'  \item{CIUpper}{The upper bound of the confidence interval for \code{Estimate}, which is derived from the normal distribution.}
#'  \item{Alpha}{The significance level used in constructing the confidence interval.}
#'  \item{Pvalue}{P-value associated with \code{Estimate}, which is derived from the normal distribution.}
#'  \item{Tau2}{The variance of the balanced horizontal pleiotropy. \code{Tau2} is calculated By using all IVs in the data before conducting the IV selection.}
#'  \item{SNPs}{The number of SNPs after IV selection.}
#'  \item{Condition}{The estimated effective sample size. It is recommended to be greater than 10 for the mdIVW estimator to achieve reliable asymptotic properties. See 'Details'.}
#'
#' @examples mr_mdIVW(gamma_hat = Bx_exp, gamma_se = Bxse_exp, Gamma_hat = By_exp, Gamma_se = Byse_exp)
#'
#' @export


mr_mdIVW <- function(gamma_hat,gamma_se,Gamma_hat,Gamma_se, rs_hat=NULL,rs_se=NULL, over.dispersion=TRUE, delta=0, alpha=0.05){
  # 用mdVIW估计水平多效性
  if(over.dispersion==TRUE){
    miu1_hat <- sum(gamma_hat*Gamma_hat/Gamma_se^2)
    miu2_hat <- sum((gamma_hat^2-gamma_se^2)/Gamma_se^2)
    var_miu2_hat <- sum((4*gamma_hat^2*gamma_se^2-2*gamma_se^4)/Gamma_se^4)
    cov_miu12_hat <- 2*sum((gamma_hat*Gamma_hat*gamma_se^2)/Gamma_se^4)
    beta_m<- miu1_hat/miu2_hat * (1 + cov_miu12_hat/(miu1_hat*miu2_hat) - var_miu2_hat/miu2_hat^2)

    mu2p <- miu2_hat/2 + sign(miu2_hat) * sqrt(miu2_hat^2/4 +  var_miu2_hat)
    mu1p <- miu1_hat + (cov_miu12_hat/var_miu2_hat)*(mu2p-miu2_hat)
    beta_pIVW = mu1p/mu2p

    # 用mdVIW估计水平多效性
    tau_square<- sum(((Gamma_hat-beta_m*gamma_hat)^2- Gamma_se^2-beta_m^2*gamma_se^2)*Gamma_se^(-2))/sum(Gamma_se^(-2))

    tau_square <- max(tau_square,0)
  } else {tau_square = 0}




  # 处理筛选情况
  if(delta!=0){
    prob_hat <- pnorm(rs_hat/rs_se-delta)+pnorm(-rs_hat/rs_se-delta)
    extra <-sum((gamma_hat^(4)*gamma_se^(-4) - 6*gamma_hat^(2)*gamma_se^(-2) + 3)*prob_hat*(1-prob_hat))
    sel = (abs(rs_hat)/rs_se) > delta
  } else {extra <- 0;sel=rep(TRUE,length(gamma_hat))}


  # 经过筛选之后的数据集
  gamma_hat <- gamma_hat[sel]
  Gamma_hat <- Gamma_hat[sel]
  gamma_se <- gamma_se[sel]
  Gamma_se <- Gamma_se[sel]

  # 有效样本量
  p_hat <- sum(sel)
  kappa_hat <- sum((gamma_hat^2-gamma_se^2)/gamma_se^2) / p_hat
  phi_hat <- sqrt(extra/p_hat)
  ### extra值可能为负，pha_hat=NA
  eta <- kappa_hat*sqrt(p_hat)/max(1,phi_hat,na.rm = TRUE)

  # 统计量 dIVW mdIVW
  miu1_hat <- sum(gamma_hat*Gamma_hat/Gamma_se^2)
  miu2_hat <- sum((gamma_hat^2-gamma_se^2)/Gamma_se^2)
  #* 协方差
  var_miu1_hat <- sum((Gamma_hat^2*gamma_se^2+gamma_hat^2*(Gamma_se^2+tau_square)-gamma_se^2*(Gamma_se^2+tau_square))/Gamma_se^4)
  var_miu2_hat <- sum((4*gamma_hat^2*gamma_se^2-2*gamma_se^4)/Gamma_se^4)
  cov_miu12_hat <- 2*sum((gamma_hat*Gamma_hat*gamma_se^2)/Gamma_se^4)

  # #*dIVW----------------------------------------------------------------------------------
  # #*点估计
  beta_d <- miu1_hat/miu2_hat
  # #*方差
  w <- gamma_hat^2/Gamma_se^2
  v <- gamma_se^2 /Gamma_se^2
  # V_dIVW <- sum(w*(1+tau_square/Gamma_se^2)+beta_d^2*v*(w+v))/sum(w-v)^2
  # se_d <- sqrt(V_dIVW)

  #*mdIVW---------------------------------------------------------------------------------
  #*点估计
  mod_md <- (1 + cov_miu12_hat/(miu1_hat*miu2_hat) - var_miu2_hat/miu2_hat^2)
  beta_m <- beta_d * mod_md
  #*方差
  a1 <- sum(gamma_se^6*Gamma_se^(-6)*(6*gamma_se^(-2)*gamma_hat^2+8))
  a2 <- sum(Gamma_se^(-4)*gamma_se^4*(gamma_se^(-2)*gamma_hat^2+1)+Gamma_se^(-6)*gamma_se^4*tau_square*(gamma_se^(-2)*gamma_hat^2+1))

  diff<- 2*miu2_hat^(-4)*(var_miu1_hat*var_miu2_hat
                          - 6*beta_m*cov_miu12_hat*var_miu2_hat
                          + 2*cov_miu12_hat^2
                          + 3*beta_m^2*var_miu2_hat^2
                          - miu2_hat*beta_m^2*a1
                          - 2*miu2_hat*a2)
  V_mdIVW_st <- sum(w*(1+tau_square/Gamma_se^2)+beta_m^2*v*(w+v))/sum(w-v)^2
  V_mdIVW <- V_mdIVW_st-diff
  se_m <- sqrt(V_mdIVW)
  if (V_mdIVW <0 ) se_m <- sqrt(V_mdIVW_st)

  pval = 2*pnorm(abs(beta_m),0,se_m,lower.tail = FALSE)
  CI_ll = beta_m+qnorm(alpha/2)*se_m
  CI_ul = beta_m-qnorm(alpha/2)*se_m

  return(new("mdIVW",
             Over.dispersion = over.dispersion,
             Delta = delta,
             Estimate = beta_m,
             StdError = se_m,
             CILower = CI_ll,
             CIUpper = CI_ul,
             Alpha = alpha,
             Pvalue = as.numeric(pval),
             Tau2 = tau_square,
             SNPs = length(Gamma_hat),
             Condition = eta)
  )

}
