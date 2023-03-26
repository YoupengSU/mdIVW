mr_mdIVW <- function(gamma_hat,Gamma_hat,gamma_se,Gamma_se,rs_hat=NULL,rs_se=NULL,lamda=0,sig.level=0.95,pleiotropy=FALSE){
  # 计算水平多效性
  if(pleiotropy==TRUE){
    miu1_hat <- sum(gamma_hat*Gamma_hat/Gamma_se^2)
    miu2_hat <- sum((gamma_hat^2-gamma_se^2)/Gamma_se^2)
    var_miu2_hat <- sum((4*gamma_hat^2*gamma_se^2-2*gamma_se^4)/Gamma_se^4)
    cov_miu12_hat <- 2*sum((gamma_hat*Gamma_hat*gamma_se^2)/Gamma_se^4)
    beta_tuta<- miu1_hat/miu2_hat * (1 + cov_miu12_hat/(miu1_hat*miu2_hat) - var_miu2_hat/miu2_hat^2)
    
    tau_square<- sum(((Gamma_hat-beta_tuta*gamma_hat)^2- 
                        Gamma_se^2-beta_tuta^2*gamma_se^2)*Gamma_se^(-2))/sum(Gamma_se^(-2))
    tau_square <- max(tau_square,0)
  } else tau_square = 0
  
  # 处理筛选情况
  if(lamda!=0){
    prob_hat <- pnorm(rs_hat/rs_se-lamda)+pnorm(-rs_hat/rs_se-lamda)
    extra <-sum((gamma_hat^(4)*gamma_se^(-4) - 6*gamma_hat^(2)*gamma_se^(-2) + 3)*prob_hat*(1-prob_hat))
    sel = (abs(rs_hat)/rs_se) > lamda
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
  eta <- kappa_hat*sqrt(p_hat)/max(1,phi_hat)
  
  # 统计量 dIVW mdIVW
  miu1_hat <- sum(gamma_hat*Gamma_hat/Gamma_se^2)
  miu2_hat <- sum((gamma_hat^2-gamma_se^2)/Gamma_se^2)
  #* 协方差
  var_miu1_hat <- sum((Gamma_hat^2*gamma_se^2+gamma_hat^2*(Gamma_se^2+tau_square)-gamma_se^2*(Gamma_se^2+tau_square))/Gamma_se^4)
  var_miu2_hat <- sum((4*gamma_hat^2*gamma_se^2-2*gamma_se^4)/Gamma_se^4)
  cov_miu12_hat <- 2*sum((gamma_hat*Gamma_hat*gamma_se^2)/Gamma_se^4)
  #* 估计量
  beta_hat <- miu1_hat/miu2_hat
  beta_tuta<- miu1_hat/miu2_hat * (1 + cov_miu12_hat/(miu1_hat*miu2_hat) - var_miu2_hat/miu2_hat^2)
  
  #* 估计方差
  a1 <- sum(gamma_se^6*Gamma_se^(-6)*(6*gamma_se^(-2)*gamma_hat^2+2))
  a2 <- sum(Gamma_se^(-4)*gamma_se^2*gamma_hat^2+Gamma_se^(-6)*gamma_se^2*tau_square*gamma_hat^2)
  
  diff<- 2*miu2_hat^(-4)*(var_miu1_hat*var_miu2_hat
                          - 6*beta_tuta*cov_miu12_hat*var_miu2_hat
                          + 2*cov_miu12_hat^2 
                          + 3*beta_tuta^2*var_miu2_hat^2
                          - miu2_hat*beta_tuta^2*a1
                          - 2*miu2_hat*a2)
  
  
  w <- gamma_hat^2/Gamma_se^2
  v <- gamma_se^2 /Gamma_se^2
  #* 一阶方差
  V_dIVW_st <- sum(w*(1+tau_square/Gamma_se^2)+beta_tuta^2*v*(w+v))/sum(w-v)^2
  #* mdIVW 方差
  Variance <- V_dIVW_st-diff
  se <- sqrt(Variance)
  d_se <- sqrt(V_dIVW_st)
  return(c(beta_hat,d_se,beta_tuta,se,tau_square,eta,kappa_hat))
}
