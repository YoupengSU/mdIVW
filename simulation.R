rm(list = ls())

# Function to caculate the dIVW, mdIVW and piVW estiamte
mr_dmp <- function(gamma_hat,Gamma_hat,gamma_se,Gamma_se,rs_hat=NULL,rs_se=NULL,delta=0,sig.level=0.95,pleiotropy=FALSE){
  
  # estimate the pleiotropy
  if(pleiotropy==TRUE){                                         
    miu1_hat <- sum(gamma_hat*Gamma_hat/Gamma_se^2)
    miu2_hat <- sum((gamma_hat^2-gamma_se^2)/Gamma_se^2)
    var_miu2_hat <- sum((4*gamma_hat^2*gamma_se^2-2*gamma_se^4)/Gamma_se^4)
    cov_miu12_hat <- 2*sum((gamma_hat*Gamma_hat*gamma_se^2)/Gamma_se^4)
    beta_m<- miu1_hat/miu2_hat * (1 + cov_miu12_hat/(miu1_hat*miu2_hat) - var_miu2_hat/miu2_hat^2)
    
    mu2p <- miu2_hat/2 + sign(miu2_hat) * sqrt(miu2_hat^2/4 +  var_miu2_hat)
    mu1p <- miu1_hat + (cov_miu12_hat/var_miu2_hat)*(mu2p-miu2_hat)
    beta_pIVW = mu1p/mu2p
    
    #* using the mdIVW estimator to estimate the pleiotropy
    tau_square<- sum(((Gamma_hat-beta_m*gamma_hat)^2- Gamma_se^2-beta_m^2*gamma_se^2)*Gamma_se^(-2))/sum(Gamma_se^(-2))
    
    #* using the pIVW estimator to estimate the pleiotropy
    tau2 = sum(((Gamma_hat - beta_pIVW * gamma_hat)^2 - Gamma_se^2 - beta_pIVW^2 * gamma_se^2) * Gamma_se^(-2))/sum(Gamma_se^(-2))
    
    tau_square <- max(tau_square,0)
    tau2 <- max(tau2,0)
  } else {tau_square = 0; tau2=0}

  
  # the case with IV selection
  if(delta!=0){
    prob_hat <- pnorm(rs_hat/rs_se-delta)+pnorm(-rs_hat/rs_se-delta)
    extra <-sum((gamma_hat^(4)*gamma_se^(-4) - 6*gamma_hat^(2)*gamma_se^(-2) + 3)*prob_hat*(1-prob_hat))
    sel = (abs(rs_hat)/rs_se) > delta
  } else {extra <- 0;sel=rep(TRUE,length(gamma_hat))}
  
  
  # data after IV selection
  gamma_hat <- gamma_hat[sel]
  Gamma_hat <- Gamma_hat[sel]
  gamma_se <- gamma_se[sel]
  Gamma_se <- Gamma_se[sel]
  
  # estimated effective sample size
  p_hat <- sum(sel)
  kappa_hat <- sum((gamma_hat^2-gamma_se^2)/gamma_se^2) / p_hat 
  phi_hat <- sqrt(extra/p_hat) 
  eta <- kappa_hat*sqrt(p_hat)/max(1,phi_hat,na.rm = TRUE)  
  

  miu1_hat <- sum(gamma_hat*Gamma_hat/Gamma_se^2)
  miu2_hat <- sum((gamma_hat^2-gamma_se^2)/Gamma_se^2)

  var_miu1_hat <- sum((Gamma_hat^2*gamma_se^2+gamma_hat^2*(Gamma_se^2+tau_square)-gamma_se^2*(Gamma_se^2+tau_square))/Gamma_se^4)
  var_miu2_hat <- sum((4*gamma_hat^2*gamma_se^2-2*gamma_se^4)/Gamma_se^4)
  cov_miu12_hat <- 2*sum((gamma_hat*Gamma_hat*gamma_se^2)/Gamma_se^4)
  
  #*dIVW----------------------------------------------------------------------------------
  #* point estimate
  beta_d <- miu1_hat/miu2_hat
  #* variance estimate
  w <- gamma_hat^2/Gamma_se^2
  v <- gamma_se^2 /Gamma_se^2
  V_dIVW <- sum(w*(1+tau_square/Gamma_se^2)+beta_d^2*v*(w+v))/sum(w-v)^2
  se_d <- sqrt(V_dIVW)
  
  #*mdIVW---------------------------------------------------------------------------------
  #* point estimate
  mod_md <- (1 + cov_miu12_hat/(miu1_hat*miu2_hat) - var_miu2_hat/miu2_hat^2)
  beta_m <- beta_d * mod_md
  #* variance estimate
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
  
  #*pIVW----------------------------------------------------------------------------------
  #* point estimate
  mu2p <- miu2_hat/2 + sign(miu2_hat) * sqrt(miu2_hat^2/4 +  var_miu2_hat)
  mu1p <- miu1_hat + (cov_miu12_hat/var_miu2_hat)*(mu2p-miu2_hat)
  beta_p = mu1p/mu2p
  #* variance estimate
  se_p = 1/abs(mu2p) * sqrt(sum((gamma_hat^2/Gamma_se^2) * (1 + tau2 * Gamma_se^(-2)) 
                                + beta_p^2 * (gamma_se^2/Gamma_se^4) * (gamma_hat^2 +  gamma_se^2)))
  #*IVW-----------------------------------------------------------------------------------
  beta_IVW <-miu1_hat/sum(w) 
  V_IVW <- sum(w*(1+tau_square/Gamma_se^2)+beta_IVW^2*v*(w+v))/sum(w)^2
  se_IVW <- sqrt(V_IVW)
 
  return(c(beta_d=beta_d,beta_m=beta_m,beta_p=beta_p,beta_IVW=beta_IVW,
           se_d=se_d,se_m=se_m,se_p=se_p,se_IVW =se_IVW,
           eta=eta,kappa=kappa_hat,theta2 = miu2_hat,mod=mod_md,
           tau_m = tau_square, tau_p = tau2))
}

# Function to generate the parameters
sim_xu <- function(J,pop,beta,nx,ny,sx,se_r,V_ex,V_ey,V_U,myseed){
  # setting seed
  set.seed(myseed)
  
  sJ= round(pop*J)        # number of weak IVs
  r_e = rnorm(sJ,0, se_r) # exposure-SNP associations for weak IVs
  r = c(r_e,rep(0,J-sJ))  # exposure-SNP associations for all IVs
  
  maf = runif(J,0.1,0.5)  # minimum allele frequency
  var_SNP = 2*maf*(1-maf) # SNP variance
  
  I = beta*r # outcome-SNP associations for all IVs
  
  var.x = sum(r^2*var_SNP) + V_U + V_ex # variance of the exposure
  var.y = sum(I^2*var_SNP) + (1+beta)^2*V_U + beta^2*V_ex + V_ey # variance of the outcome
  
  r_se = sqrt((var.x - r^2*var_SNP)/(nx*var_SNP)) # standard errors for exposure-SNP associations
  I_se = sqrt((var.y - I^2*var_SNP)/(ny*var_SNP)) # standard errors for outcome-SNP associations
  
  # selection GWAS data
  rs_se.pop = sqrt(nx/sx) 
  rs_se = r_se*rs_se.pop  
  
  # combine the true values
  data_true <- as.data.frame(cbind(r,I,r_se,I_se,rs_se)) 
  return(data_true)
}


# Simulation results for once sampling
sim_parallel <- function(tau,delta,data_true){
  flag <- ifelse(tau==0,FALSE,TRUE)
  
  a = rnorm(1000,0,tau)
  Gamma = data_true[,2]+a
  
  # generate exposure/outcome-SNP association estimates
  gamma_hat = rnorm(1000,data_true[,1], data_true[,3])
  Gamma_hat = rnorm(1000,Gamma, data_true[,4])
  rs_hat = rnorm(1000,data_true[,1], data_true[,5])
  
  # generate standard errors
  gamma_se = data_true[,3]
  Gamma_se = data_true[,4]
  rs_se <- data_true[,5]
  
  # number of IVs after selection
  nsel_1 = sum((abs(rs_hat)/data_true[,5]) > delta[2])
  
  
  sel.pval = 2*pnorm( (abs(rs_hat)/rs_se),lower.tail = FALSE)
  
  # result without IV selection
  out_mdIVW_0 <- mr_dmp(gamma_hat,Gamma_hat,gamma_se,Gamma_se,rs_hat,rs_se,0,0.95,flag)
  
  # result with IV selection
  if (nsel_1>2){
    out_mdIVW_1 <- mr_dmp(gamma_hat,Gamma_hat,gamma_se,Gamma_se,rs_hat,rs_se,delta[2],0.95,flag)
  } else {out_mdIVW_1 <- rep(NA,14)}
  
  
  out <- cbind(rbind(out_mdIVW_0,out_mdIVW_1),
               delta)  
  
  return(out)
}

# Simulation results of each parameter combination based on 10,000 repetitions
summary_out <- function(data_sample,beta,Ind){
  obs <- nrow(data_sample) 
  if (Ind == 99) {data_sample <- data_sample[(floor(obs*0.005)+1):floor(obs*0.995),]
  }else if (Ind == 95){data_sample <- data_sample[(floor(obs*0.025)+1):floor(obs*0.975),]
  }else if (Ind == 90) {data_sample <- data_sample[(floor(obs*0.05)+1):floor(obs*0.95),]}
  
  data_sample <- na.omit(data_sample)
  n_obs <- nrow(data_sample)
  z <- qnorm(0.025, lower.tail = FALSE)
  # mean
  mean_d <- mean(data_sample$beta_d)
  mean_m <- mean(data_sample$beta_m)
  mean_p <- mean(data_sample$beta_p)
  # relative bias
  bias_d <- (mean_d-beta)/beta*100
  bias_m <- (mean_m-beta)/beta*100
  bias_p <- (mean_p-beta)/beta*100
  # standard deviation of 10,000 point estimates
  SD_d <- sqrt(var(data_sample$beta_d))
  SD_m <- sqrt(var(data_sample$beta_m))
  SD_p <- sqrt(var(data_sample$beta_p))
  # mean of standard errors in 10,000 repetitions
  SE_d <- mean(data_sample$se_d)
  SE_m <- mean(data_sample$se_m)
  SE_p <- mean(data_sample$se_p)

  # MSE
  MSE_d <- mean((data_sample$beta_d-beta)^2,na.rm = TRUE)
  MSE_m <- mean((data_sample$beta_m-beta)^2,na.rm = TRUE)
  MSE_p <- mean((data_sample$beta_p-beta)^2,na.rm = TRUE)

  # CP
  CP_d <- sum(data_sample$beta_d-z*data_sample$se_d<beta & data_sample$beta_d+z*data_sample$se_d > beta)/n_obs
  CP_m <- sum(data_sample$beta_m-z*data_sample$se_m<beta & data_sample$beta_m+z*data_sample$se_m > beta)/n_obs
  CP_p <- sum(data_sample$beta_p-z*data_sample$se_p<beta & data_sample$beta_p+z*data_sample$se_p > beta)/n_obs
  
  # mean effective sample sizes in 10,000 repetitions
  eta <- mean(data_sample$eta)
  kappa <- mean(data_sample$kappa)
  tau_mean_m <- mean(data_sample$tau_m)
  tau_mean_p <- mean(data_sample$tau_p)
  
  return(c(mean_d=mean_d, mean_m=mean_m, mean_p=mean_p,
           bias_d=bias_d, bias_m=bias_m, bias_p=bias_p,
           SD_d=SD_d, SD_m=SD_m, SD_p=SD_p,
           SE_d=SE_d, SE_m=SE_m, SE_p=SE_p,
           MSE_d= MSE_d, MSE_m= MSE_m, MSE_p= MSE_p,
           CP_d=CP_d, CP_m=CP_m, CP_P=CP_p,
           eta=eta, kappa=kappa,
           tau_mean_m = tau_mean_m, tau_mean_p = tau_mean_p))
}


#----------------------------Here is a test for the function defined above----------------

data_true <- sim_xu(J=1000,pop=0.05,beta=0.5,nx=150000,ny=225000,sx=75000,se_r=sqrt(5)/100,
                    V_ex = 2, V_ey = 2, V_U = 2,myseed = 1)

s <- Sys.time()
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
out <- foreach(x=1:10000,.combine='rbind') %dopar% 
  sim_parallel(tau=0.01,delta=c(0,3.7),data_true = data_true)
stopCluster(cl)
e <- Sys.time()
e-s

out <- as.data.frame(out)
out1 <- out[out$delta==0,]
out1_sort <- out1[order(out1$eta),]

kappa_0 <- sum(data_true$r^2/data_true$r_se^2)/nrow(data_true)
eta_0 <- kappa_0*sqrt(nrow(data_true))

hist(out1_sort$eta,breaks = 30)

summary(out1_sort$eta)


out2 <- out[out$delta==3.7,]
out2 <- na.omit(out2)
out2_sort <- out2[order(out2$eta),]

result <- rbind(summary_out(out1_sort,0.5,Ind=100),
                summary_out(out1_sort,0.5,Ind=99),
                summary_out(out1_sort,0.5,Ind=95),
                summary_out(out1_sort,0.5,Ind=90),
                
                summary_out(out2_sort,0.5,Ind=100),
                summary_out(out2_sort,0.5,Ind=99),
                summary_out(out2_sort,0.5,Ind=95),
                summary_out(out2_sort,0.5,Ind=90))




result <- cbind(result,delta=rep(c(0,3.7),each=4),trim=rep(c("t100","t99","t95","t90"),2))

# true values of eta and kappa
#* no IV selection
kappa_0 <- sum(data_true$r^2/data_true$r_se^2)/nrow(data_true)
eta_0 <- kappa_0*sqrt(nrow(data_true))

#* IV selection
q_delta <- pnorm(data_true$r/data_true$rs_se-3.7)+pnorm(-data_true$r/data_true$rs_se-3.7)
p_delta <- sum(q_delta)

kappa_delta <- sum(data_true$r_se^(-2)*data_true$r^2*q_delta)/p_delta
omga_delta <- sqrt(sum(data_true$I_se^(-4)*data_true$r^4*q_delta*(1-q_delta))/p_delta)
eta_delta <- kappa_delta*sqrt(p_delta)/max(1,omga_delta) 

rep(c(eta_0,eta_delta),each=4)
rep(c(kappa_0,kappa_delta),each=4)

result <- cbind(result,
                eta_true=rep(c(eta_0,eta_delta),each=4),
                kappa_true=rep(c(kappa_0,kappa_delta),each=4))

#----------------------------Simulation---------------------------------------------------
#* parameter settings---------------------------------------------------------------------
IV = c(0:100)/100           # proportion of weak IVs
se_ser = sqrt(c(1:10))/100  # SD of the distribution generating exposure-SNP associations
NX = c(10:20)*10000         # the exposure GWAS sample size
par <- expand.grid(IV,se_ser,NX) # parameter combinations
pY = 0.5 # ratio of the outcome GWAS sample size to the exposure GWAS sample size
beta = 0.5 # true causal effect
myseed = 1 # seed
rep = 10000 # repetitions under each parameter combination
tau = 0 # balanced horizontal pleiotropy
# tau = 0.01 #* the case with balanced horizontal pleiotropy

#* conduct simulation---------------------------------------------------------------------
library(doParallel)

set_num <- nrow(par)
out_total <- NULL

for (i in 1:set_num) {
  pop <- par[i,1]
  se_r <- par[i,2]
  NX <- par[i,3]
  NY <- NX*pY
  data <- sim_xu(J=1000,pop=pop,beta=beta,
                 nx=NX,ny=NY,sx=NY,se_r=se_r,
                 V_ex = 2, V_ey = 2, V_U = 2, 
                 myseed = myseed)
  
  cl <- makeCluster(8)
  registerDoParallel(cl)
  out <- foreach(x=1:rep,.combine='rbind') %dopar% 
    sim_parallel(tau=tau,delta=c(0,3.7),data_true = data)
  stopCluster(cl)
  
  out <- as.data.frame(out)
  out1 <- out[out$delta==0,]
  out1_sort <- out1[order(out1$eta),]
  
  out2 <- out[out$delta==3.7,]
  out2 <- na.omit(out2)
  out2_sort <- out2[order(out2$eta),]
  
  result <- rbind(summary_out(out1_sort,0.5,Ind=100),
                  summary_out(out1_sort,0.5,Ind=99),
                  
                  summary_out(out2_sort,0.5,Ind=100),
                  summary_out(out2_sort,0.5,Ind=99))
  
  
  result <- cbind(result,delta=rep(c(0,3.7),each=2),trim=rep(c("t100","t99"),2))
  
  # true values of eta and kappa
  #* no IV selection
  kappa_0 <- sum(data$r^2/data$r_se^2)/nrow(data)
  eta_0 <- kappa_0*sqrt(nrow(data))
  
  #* IV selection
  q_delta <- pnorm(data$r/data$rs_se-3.7)+pnorm(-data$r/data$rs_se-3.7)
  p_delta <- sum(q_delta)
  
  kappa_delta <- sum(data$r_se^(-2)*data$r^2*q_delta)/p_delta
  omga_delta <- sqrt(sum(data$I_se^(-4)*data$r^4*q_delta*(1-q_delta))/p_delta)
  eta_delta <- kappa_delta*sqrt(p_delta)/max(1,omga_delta) 
  
  result <- cbind(result,
                  eta_true=rep(c(eta_0,eta_delta),each=2),
                  kappa_true=rep(c(kappa_0,kappa_delta),each=2))
  
  out_total <- rbind(out_total,result)
  if (i%%20 ==0) cat(i,"\n")
}

# write the simulation results to CSV files
# result with tau = 0
write.csv(out_total,file = "out_total_tau0.csv")
# result wiht tau=0.01
write.csv(out_total,file = "out_total_tau0.csv")

#----------------------------PLOT---------------------------------------------------------
rm(list=ls())
# setwd("./beta05")
getwd()
# load the simulation result with no pleiotropy
data <- read.csv(file = "out_total_tau0.csv")

# load the simulation result with balanced horizontal pleiotropy
#data <- read.csv(file = "out_total_tau1.csv")

library(ggplot2)


#----------------------------*the scenario without IV selection---------------------------
data1 <- data[data$delta==0,]
data_ana <-  data1[data1$trim=="t99",]

#**bias------------------------------------------------------------------------------------

data_bias <- data_ana[c("bias_d","bias_m","bias_p","eta","kappa","eta_true","kappa_true")]
data_w <- data_bias[data_bias$eta>5,]


bias_d <- data.frame(x=data_w$eta,y=data_w$bias_d,cat="dIVW")
bias_m <- data.frame(x=data_w$eta,y=data_w$bias_m,cat="mdIVW")
bias_p <- data.frame(x=data_w$eta,y=data_w$bias_p,cat="pIVW")

# mdIVW VS dIVW VS pIVW
bias <- rbind(bias_d,bias_p,bias_m)
bias$cat <- factor(bias$cat,levels=c("dIVW","pIVW","mdIVW"))


ggplot(bias,aes(x,y,color=cat))+
  geom_point(size=0.8)+
  scale_x_continuous(limits = c(5,30),breaks = c(10,20,30))+
  # scale_y_continuous(limits = c(-10,13), breaks = c(-10,-5,0,5,10,15,20))+
  labs(x="",y="",color="") +
  scale_color_manual(values = c("dIVW"="red","pIVW"="orange","mdIVW"="blue")) +
  # geom_hline(yintercept = -1, colour="blue",linetype=2)+
  # geom_hline(yintercept = 1, colour="blue",linetype=2)+
  geom_hline(yintercept = 0, colour="blue",linetype=2)+
  #theme_bw()+
  theme(legend.position = c(.85,.85),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(0.45,"cm"),
        legend.key.height = unit(0.4,"cm")) +
  guides(color=guide_legend(override.aes = list(size = 2.5))) 

#**SD----------------------------------------------------------------------------------
data_ana <-  data1[data1$trim=="t99",]
data_SD <- data_ana[c("SD_d","SD_m","SD_p","SE_d","SE_m","SE_p","eta","kappa","eta_true","kappa_true")]
data_w <- data_SD[data_SD$eta>5.05,]

SD_d <- data.frame(x=data_w$eta,y=data_w$SD_d,cat="dIVW")
SD_m <- data.frame(x=data_w$eta,y=data_w$SD_m,cat="mdIVW")
SD_p <- data.frame(x=data_w$eta,y=data_w$SD_p,cat="pIVW")

# mdIVW VS dIVW VS pIVW
SD <- rbind(SD_d,SD_p,SD_m)
SD$cat <- factor(SD$cat, levels = c("dIVW","pIVW","mdIVW"))
ggplot(SD,aes(x,y,color=cat))+
  geom_point(size=0.8)+
  scale_x_continuous(limits = c(5,30),breaks = c(10,20,30))+
  # scale_y_continuous(limits = c(0.08,0.7), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6))+
  labs(x="",y="",color="") +
  scale_color_manual(values = c("dIVW"="red","pIVW"="orange","mdIVW"="blue")) +
  #theme_bw()+
  theme(legend.position = c(.85,.85),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(0.45,"cm"),
        legend.key.height = unit(0.4,"cm")) +
  guides(color=guide_legend(override.aes = list(size = 2.5)))

#**MSE--------------------------------------------------------------------------------------
data_MSE <- data_ana[c("MSE_d","MSE_m","MSE_p","eta","kappa","eta_true","kappa_true")]
data_w <- data_MSE[data_MSE$eta>5.05,]

MSE_d <- data.frame(x=data_w$eta,y=data_w$MSE_d,cat="dIVW")
MSE_m <- data.frame(x=data_w$eta,y=data_w$MSE_m,cat="mdIVW")
MSE_p <- data.frame(x=data_w$eta,y=data_w$MSE_p,cat="pIVW")

# mdIVW VS dIVW VS pIVW
MSE <- rbind(MSE_d,MSE_p,MSE_m)
MSE$cat <- factor(MSE$cat, levels = c("dIVW","pIVW","mdIVW"))
ggplot(MSE,aes(x,y,color=cat))+
  geom_point(size=0.8)+
  scale_x_continuous(limits = c(5,30),breaks = c(10,20,30))+
  #scale_y_continuous(limits = c(0,0.22))+
  labs(x="",y="",color="") +
  scale_color_manual(values = c("dIVW"="red","pIVW"="orange","mdIVW"="blue")) +
  #theme_bw()+
  theme(legend.position = c(.85,.85),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(0.45,"cm"),
        legend.key.height = unit(0.4,"cm")) +
  guides(color=guide_legend(override.aes = list(size = 2.5)))

#**CP--------------------------------------------------------------------------------------

data_ana <- data1[data1$trim=="t99",]

data_CP <- data_ana[c("CP_d","CP_m","CP_P","eta","kappa","eta_true","kappa_true")]
data_w <- data_CP[data_CP$eta>5,]

CP_d <- data.frame(x=data_w$eta,y=data_w$CP_d,cat="dIVW")
CP_m <- data.frame(x=data_w$eta,y=data_w$CP_m,cat="mdIVW")
CP_p <- data.frame(x=data_w$eta,y=data_w$CP_P,cat="pIVW")

CP <- rbind(CP_d,CP_p,CP_m)
CP$cat <- factor(CP$cat,levels = c("dIVW","pIVW","mdIVW"))
ggplot(CP,aes(x,y,color=cat))+
  geom_point(size=0.8)+
  scale_x_continuous(limits = c(5,30),breaks = c(10,20,30))+
  scale_y_continuous(limits = c(0.90,1),breaks = c(0.90,0.95,1.00))+
  labs(x="",y="",color="") +
  # scale_color_manual(values = c("dIVW"="red","pIVW"="#013565","mdIVW"="blue")) +
  scale_color_manual(values = c("dIVW"="red","pIVW"="orange","mdIVW"="blue")) +
  #scale_color_discrete(name="cat",labels=c("dIVW","mdIVW"))+
  geom_hline(yintercept = 0.95, colour="blue",linetype=2)+
  #theme_bw()+
  theme(legend.position = c(.85,.85),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(0.45,"cm"),
        legend.key.height = unit(0.4,"cm")) +
  guides(color=guide_legend(override.aes = list(size = 2.5)))


#----------------------------*the scenario with IV selection------------------------------
data2 <- data[data$delta==3.7,]
data_ana <-  data2[data2$trim=="t99",]

#**bias------------------------------------------------------------------------------------

data_bias <- data_ana[c("bias_d","bias_m","bias_p","eta","kappa","eta_true","kappa_true")]
data_w <- data_bias[data_bias$eta>5,]


bias_d <- data.frame(x=data_w$eta,y=data_w$bias_d,cat="dIVW")
bias_m <- data.frame(x=data_w$eta,y=data_w$bias_m,cat="mdIVW")
bias_p <- data.frame(x=data_w$eta,y=data_w$bias_p,cat="pIVW")

# mdIVW VS dIVW VS pIVW
bias <- rbind(bias_d,bias_p,bias_m)
bias$cat <- factor(bias$cat,levels=c("dIVW","pIVW","mdIVW"))


ggplot(bias,aes(x,y,color=cat))+
  geom_point(size=0.8)+
  scale_x_continuous(limits = c(5,30),breaks = c(10,20,30))+
  # scale_y_continuous(limits = c(-1,2), breaks = c(-1,0,1,2))+
  labs(x="",y="",color="") +
  scale_color_manual(values = c("dIVW"="red","pIVW"="orange","mdIVW"="blue")) +
  geom_hline(yintercept = 0, colour="blue",linetype=2)+
  #theme_bw()+
  theme(legend.position = c(.85,.85),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(0.45,"cm"),
        legend.key.height = unit(0.4,"cm")) +
  guides(color=guide_legend(override.aes = list(size = 2.5))) 

#**SD-------------------------------------------------------------------------------------
data_SD <- data_ana[c("SD_d","SD_m","SD_p","SE_d","SE_m","SE_p","eta","kappa","eta_true","kappa_true")]
data_w <- data_SD[data_SD$eta>5,]

SD_d <- data.frame(x=data_w$eta,y=data_w$SD_d,cat="dIVW")
SD_m <- data.frame(x=data_w$eta,y=data_w$SD_m,cat="mdIVW")
SD_p <- data.frame(x=data_w$eta,y=data_w$SD_p,cat="pIVW")

# mdIVW VS dIVW VS pIVW
SD <- rbind(SD_d,SD_p,SD_m)
SD$cat <- factor(SD$cat, levels = c("dIVW","pIVW","mdIVW"))
ggplot(SD,aes(x,y,color=cat))+
  geom_point(size=0.8)+
  scale_x_continuous(limits = c(4,30),breaks = c(10,20,30))+
  # scale_y_continuous(limits = c(0.02,0.3), breaks = c(0,0.1,0.2,0.3))+
  labs(x="",y="",color="") +
  scale_color_manual(values = c("dIVW"="red","pIVW"="orange","mdIVW"="blue")) +
  #theme_bw()+
  theme(legend.position = c(.85,.85),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(0.45,"cm"),
        legend.key.height = unit(0.4,"cm")) +
  guides(color=guide_legend(override.aes = list(size = 2.5)))

#**MSE--------------------------------------------------------------------------------------
data_MSE <- data_ana[c("MSE_d","MSE_m","MSE_p","eta","kappa","eta_true","kappa_true")]
data_w <- data_MSE[data_MSE$eta>5,]

MSE_d <- data.frame(x=data_w$eta,y=data_w$MSE_d,cat="dIVW")
MSE_m <- data.frame(x=data_w$eta,y=data_w$MSE_m,cat="mdIVW")
MSE_p <- data.frame(x=data_w$eta,y=data_w$MSE_p,cat="pIVW")

# mdIVW VS dIVW VS pIVW
MSE <- rbind(MSE_d,MSE_p,MSE_m)
MSE$cat <- factor(MSE$cat, levels = c("dIVW","pIVW","mdIVW"))
ggplot(MSE,aes(x,y,color=cat))+
  geom_point(size=0.8)+
  scale_x_continuous(limits = c(4,30),breaks = c(10,20,30))+
  # scale_y_continuous(limits = c(0,0.22))+
  labs(x="",y="",color="") +
  scale_color_manual(values = c("dIVW"="red","pIVW"="orange","mdIVW"="blue")) +
  #theme_bw()+
  theme(legend.position = c(.85,.85),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(0.45,"cm"),
        legend.key.height = unit(0.4,"cm")) +
  guides(color=guide_legend(override.aes = list(size = 2.5)))

#**CP--------------------------------------------------------------------------------------

data_ana <- data2[data2$trim=="t100",]

data_CP <- data_ana[c("CP_d","CP_m","CP_P","eta","kappa","eta_true","kappa_true")]
data_w <- data_CP[data_CP$eta>5,]

CP_d <- data.frame(x=data_w$eta,y=data_w$CP_d,cat="dIVW")
CP_m <- data.frame(x=data_w$eta,y=data_w$CP_m,cat="mdIVW")
CP_p <- data.frame(x=data_w$eta,y=data_w$CP_P,cat="pIVW")
CP <- rbind(CP_d,CP_p,CP_m)
CP$cat <- factor(CP$cat,levels = c("dIVW","pIVW","mdIVW"))
ggplot(CP,aes(x,y,color=cat))+
  geom_point(size=0.8)+
  scale_x_continuous(limits = c(4,30),breaks = c(10,20,30))+
  scale_y_continuous(limits = c(0.90,1),breaks = c(0.90,0.95,1.00))+
  labs(x="",y="",color="") +
  # scale_color_manual(values = c("dIVW"="red","pIVW"="#013565","mdIVW"="blue")) +
  scale_color_manual(values = c("dIVW"="red","pIVW"="orange","mdIVW"="blue")) +
  #scale_color_discrete(name="cat",labels=c("dIVW","mdIVW"))+
  geom_hline(yintercept = 0.95, colour="blue",linetype=2)+
  #theme_bw()+
  theme(legend.position = c(.85,.85),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        legend.key.width = unit(0.45,"cm"),
        legend.key.height = unit(0.4,"cm")) +
  guides(color=guide_legend(override.aes = list(size = 2.5)))


#----------------------------Compare with outher method-----------------------------------
# tau=0.01
# delta = 3.7

# Simulation results for once sampling
sim_parallel_table <- function(tau,delta,data_true){
  library(mr.raps)
  library(MendelianRandomization)
  
  flag <- ifelse(tau==0,FALSE,TRUE)
  
  # geneate sample data
  a = rnorm(1000,0,tau)
  Gamma = data_true[,2]+a
  
  gamma_hat = rnorm(1000,data_true[,1], data_true[,3])
  Gamma_hat = rnorm(1000,Gamma, data_true[,4])
  rs_hat = rnorm(1000,data_true[,1], data_true[,5])
  
  gamma_se = data_true[,3]
  Gamma_se = data_true[,4]
  rs_se <- data_true[,5]
  
  
  
  
  #* without IV selection------------------------------------------------------------------
  if (delta==0){
    
    MRinput <- mr_input(bx = gamma_hat,
                        bxse = gamma_se,
                        by = Gamma_hat,
                        byse = Gamma_se)
    #* dIVW, mdIVW, pIVW, and IVW
    out_dmp <-  mr_dmp(gamma_hat,Gamma_hat,gamma_se,Gamma_se,rs_hat,rs_se,0,0.95,flag)
    
    #* MR-Median
    median <- mr_median(MRinput, weighting = "weighted", distribution = "normal",
                        alpha = 0.05, iterations = 1000, seed = 314159265)
    out_median <- c(median$Estimate,median$StdError)
    
    #* MR-Egger
    egger <- mr_egger(MRinput)
    out_egger <- c(egger$Estimate,egger$StdError.Est)
    
    #* MR-RAPS
    if (flag) RAPS <- {mr.raps.overdispersed.robust(b_exp=gamma_hat, 
                                                    b_out=Gamma_hat, 
                                                    se_exp=gamma_se, 
                                                    se_out=Gamma_se)
    } else {
      RAPS <- mr.raps.simple.robust(b_exp=gamma_hat, 
                                    b_out=Gamma_hat, 
                                    se_exp=gamma_se, 
                                    se_out=Gamma_se)
    }
    out_RAPS <- c(RAPS$beta.hat, RAPS$beta.se)
    
  }
  
  
  
  # with IV selection
  #* number of selected IVs
  nsel_1 = sum((abs(rs_hat)/data_true[,5]) > delta)
  
  if (delta>0 & nsel_1>2){
    
    
    ind <- abs(rs_hat)/rs_se > delta
    MRinput <- mr_input(bx = gamma_hat[ind],
                        bxse = gamma_se[ind],
                        by = Gamma_hat[ind],
                        byse = Gamma_se[ind])
    
    #* dIVW, mdIVW, and pIVW
    out_dmp <- mr_dmp(gamma_hat,Gamma_hat,gamma_se,Gamma_se,rs_hat,rs_se,delta,0.95,flag)
    
    #* MR-Median
    median <- mr_median(MRinput, weighting = "weighted", distribution = "normal",
                        alpha = 0.05, iterations = 1000, seed = 314159265)
    out_median <- c(median$Estimate,median$StdError)
    
    #* MR-Egger
    egger <- mr_egger(MRinput)
    out_egger <- c(egger$Estimate,egger$StdError.Est)
    
    
    #* MR-RAPS
    if (flag) RAPS <- {mr.raps.overdispersed.robust(b_exp=gamma_hat[ind], 
                                                    b_out=Gamma_hat[ind], 
                                                    se_exp=gamma_se[ind], 
                                                    se_out=Gamma_se[ind])
    } else  {
      RAPS <- mr.raps.simple.robust(b_exp=gamma_hat[ind], 
                                    b_out=Gamma_hat[ind], 
                                    se_exp=gamma_se[ind], 
                                    se_out=Gamma_se[ind])
    }
    
    out_RAPS <- c(RAPS$beta.hat, RAPS$beta.se)
    
    
  } else if(delta>0 & nsel_1<=2) {
    out_dmp <- rep(NA,14)
    out_median <- rep(NA,2)
    out_egger <- rep(NA,2)
    out_RAPS <- rep(NA,2)}
  
  
  
  out <- c(out_dmp,out_median,out_egger,out_RAPS)  
  return(out)
}
# sample = out_sort
# beta = 0.5
summary_out_table <- function(sample,beta){
  obs <- nrow(sample) 
  names(sample)[15:20] <- c("beta_Median","se_Median","beta_Egger","se_Egger","beta_RAPS","se_RAPS")
  # T100----------------------------------------------------------------------------------
  data_sample <- na.omit(sample)
  n_obs <- nrow(data_sample)
  z <- qnorm(0.025, lower.tail = FALSE)
  #平均值
  mean_d <- mean(data_sample$beta_d)
  mean_m <- mean(data_sample$beta_m)
  mean_p <- mean(data_sample$beta_p)
  #相对偏倚
  bias_d <- (mean_d-beta)/beta*100
  bias_m <- (mean_m-beta)/beta*100
  bias_p <- (mean_p-beta)/beta*100
  #标准差
  SD_d <- sqrt(var(data_sample$beta_d))
  SD_m <- sqrt(var(data_sample$beta_m))
  SD_p <- sqrt(var(data_sample$beta_p))
  #标准误均值
  SE_d <- mean(data_sample$se_d)
  SE_m <- mean(data_sample$se_m)
  SE_p <- mean(data_sample$se_p)
  
  #MSE均方误差
  MSE_d <- mean((data_sample$beta_d-beta)^2,na.rm = TRUE)
  MSE_m <- mean((data_sample$beta_m-beta)^2,na.rm = TRUE)
  MSE_p <- mean((data_sample$beta_p-beta)^2,na.rm = TRUE)
  
  #覆盖率
  CP_d <- sum(data_sample$beta_d-z*data_sample$se_d<beta & data_sample$beta_d+z*data_sample$se_d > beta)/n_obs
  CP_m <- sum(data_sample$beta_m-z*data_sample$se_m<beta & data_sample$beta_m+z*data_sample$se_m > beta)/n_obs
  CP_p <- sum(data_sample$beta_p-z*data_sample$se_p<beta & data_sample$beta_p+z*data_sample$se_p > beta)/n_obs
  
  #有效样本量
  eta <- mean(data_sample$eta)
  kappa <- mean(data_sample$kappa)
  
  char_dmpt100 <- matrix(c(mean_d=mean_d, mean_m=mean_m, mean_p=mean_p,
                           bias_d=bias_d, bias_m=bias_m, bias_p=bias_p,
                           SD_d=SD_d, SD_m=SD_m, SD_p=SD_p,
                           SE_d=SE_d, SE_m=SE_m, SE_p=SE_p,
                           MSE_d=MSE_d, MSE_m=MSE_m, MSE_p=MSE_p,
                           CP_d=CP_d, CP_m=CP_m, CP_P=CP_p), nrow = 3, byrow = F)
  
  
  # IVW Egger Median RAPS------------------------------------------------------------------------
  data_other <- sample[,c(4,8,15:20)]
  n_obs <- nrow(data_other)
  #平均值
  mean_IVW   <- mean(data_other$beta_IVW)
  mean_Egger <- mean(data_other$beta_Egger)
  mean_Median <- mean(sample$beta_Median)
  mean_RAPS  <- mean(data_other$beta_RAPS)
  #相对偏倚
  bias_IVW   <- (mean_IVW-beta)/beta*100
  bias_Egger <- (mean_Egger-beta)/beta*100
  bias_Median <- (mean_Median-beta)/beta*100
  bias_RAPS  <- (mean_RAPS-beta)/beta*100
  #标准差
  SD_IVW   <- sqrt(var(data_other$beta_IVW))
  SD_Egger <- sqrt(var(data_other$beta_Egger))
  SD_Median <- sqrt(var(data_other$beta_Median))
  SD_RAPS  <- sqrt(var(data_other$beta_RAPS))
  #标准误均值
  SE_IVW   <- mean(data_other$se_IVW)
  SE_Egger <- mean(data_other$se_Egger)
  SE_Median <- mean(data_other$se_Median)
  SE_RAPS  <- mean(data_other$se_RAPS)
  
  #MSE均方误差
  MSE_IVW <- mean((data_other$beta_IVW-beta)^2,na.rm = TRUE)
  MSE_Egger <- mean((data_other$beta_Egger-beta)^2,na.rm = TRUE)
  MSE_Median <- mean((data_other$beta_Median-beta)^2,na.rm = TRUE)
  MSE_RAPS <- mean((data_other$beta_RAPS-beta)^2,na.rm = TRUE)
  
  #覆盖率
  CP_IVW   <- sum(data_other$beta_IVW-z*data_other$se_IVW<beta & data_other$beta_IVW+z*data_other$se_IVW > beta)/n_obs
  CP_Egger <- sum(data_other$beta_Egger-z*data_other$se_Egger<beta & data_other$beta_Egger+z*data_other$se_Egger > beta)/n_obs
  CP_Median <- sum(data_other$beta_Median-z*data_other$se_Median<beta & data_other$beta_Median+z*data_other$se_Median > beta)/n_obs
  CP_RAPS  <- sum(data_other$beta_RAPS-z*data_other$se_RAPS<beta & data_other$beta_RAPS+z*data_other$se_RAPS > beta)/n_obs
  
  char_IER <- matrix(c(mean_IVW,mean_Median,mean_Egger,mean_RAPS,
                       bias_IVW,bias_Median,bias_Egger,bias_RAPS,
                       SD_IVW, SD_Median, SD_Egger,  SD_RAPS,
                       SE_IVW, SE_Median, SE_Egger,  SE_RAPS,
                       MSE_IVW,MSE_Median, MSE_Egger, MSE_RAPS,
                       CP_IVW, CP_Median, CP_Egger,  CP_RAPS), nrow = 4, byrow = F)
  char_IER <- rbind(char_IER,char_dmpt100)
  method <- c("IVW","MR-Median","MR-Egger","MR-RAPS","dIVW","mdIVW","pIVW")
  result <- cbind(method,char_IER,eta,kappa)
  result <- as.data.frame(result)
  names(result) <- c("method","mean","bias","SD","SE","MSE","CP","eta","kappa")
  
  return(result)
}

#----------------------------*Conduct simulation------------------------------------------
#* without IV selection and no horizontal pleiotropy--------------------------------------
data_true <- sim_xu(J=1000,pop=0.05,beta=0.5,nx=150000,ny=75000,sx=75000,
                    se_r=sqrt(5)/100,V_ex = 2, V_ey = 2, V_U = 2,myseed = 1)

time_start <- Sys.time()
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
out <- foreach(x=1:5000,.combine='rbind') %dopar% 
  sim_parallel_table(tau=0,delta=0,data_true = data_true)
stopCluster(cl)
time_end <- Sys.time()
time_end - time_start

out <- as.data.frame(out)
out <- na.omit(out)
out_sort <- out[order(out$eta),]

summary_out_table(sample = out_sort,0.5)


#* with IV selection and balanced horizontal pleiotropy-----------------------------------
data_true <- sim_xu(J=1000,pop=0.15,beta=0.5,nx=150000,ny=75000,sx=75000,
                    se_r=sqrt(5)/100,V_ex = 2, V_ey = 2, V_U = 2,myseed = 9)

time_start <- Sys.time()
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
out <- foreach(x=1:2000,.combine='rbind') %dopar% 
  sim_parallel_table(tau=0.01,delta=3.7,data_true = data_true)
stopCluster(cl)
time_end <- Sys.time()
time_end - time_start

out <- as.data.frame(out)
out <- na.omit(out)
out_sort1 <- out[order(out$eta),]

summary_out_table(sample = out_sort1,0.5)
