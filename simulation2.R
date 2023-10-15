rm(list = ls())
set.seed(123)
library(tidyverse)
library(lme4)
library(MASS)
library(foreach)
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

n_sim <- 10000
m <- 20
n_i_mean <- 10
n <- n_i_mean + 8
pi <- 0.5 # 0.5 for scenario 1 and 0.6 for scenario 2
p <- 2

tictoc::tic()
sim_result <- foreach(k = 1:n_sim, .combine = cbind, .packages = c("tidyverse", "lme4", "MASS")) %dopar% {
  S <- rbinom(m, 1, prob = 0.6)
  A <- rep(0, m)
  if(sum(S)*pi == round(sum(S)*pi)){
    A[S == 1] <- sample(c(rep(1, sum(S==1) * pi), rep(0, sum(S==1) * pi)), size = sum(S == 1), replace = F)
    A[S == 0] <- sample(c(rep(1, sum(S==0) * pi), rep(0, sum(S==0) * pi)), size = sum(S == 0), replace = F)
  } else {
    A[S == 1] <- sample(c(rep(1, sum(S==1) * pi), rep(0, sum(S==1) * pi), rbernoulli(1)), size = sum(S == 1), replace = F)
    A[S == 0] <- sample(c(rep(1, sum(S==0) * pi), rep(0, sum(S==0) * pi), rbernoulli(1)), size = sum(S == 0), replace = F)
    
  }
  sim_data <- map_dfr(1:m, function(i){
    x <- rnorm(n, sd = 2)
    y <- 2*(S[i]-0.6)*A[i] + x + 2 * mean(x) + rnorm(n, sd = 5) #+ rgamma(1, 25)
    n_i <- n
    # y <- 5 * (S[i]+x - 0.6)*A[i] +  rnorm(n, sd = 5) + rgamma(1, 5)
    # n_i <- round(runif(1, min = n_i_mean-8, max = n_i_mean+8))
    index <- sample(1:n, size = n_i, replace = F)
    data.frame(outcome = y,
               cluster = i,
               x = x,
               S = S[i],
               treatment= A[i])[index,]
  })
  cluster_level_data <- sim_data %>% group_by(cluster) %>% 
    summarise(outcome = mean(outcome),
              x = mean(x),
              S = mean(S),
              trt = mean(treatment))
  mu_x_hat <- mean(sim_data$x)
  mu_S_hat <- mean(cluster_level_data$S)
  
  tryCatch({
  
  # mix_ancova
  mix_ancova <- lmer(outcome~treatment+x+S + (1|cluster), data = sim_data, REML = F)
  zz <- summary(mix_ancova)
  beta <- zz$coefficients[,1]
  tau2 <- zz$varcor$cluster[1]
  sigma2 <- zz$sigma^2
  cl_res <- data.frame(res = sim_data$outcome - model.matrix(mix_ancova) %*% beta, 
                       cluster = sim_data$cluster, 
                       a = sim_data$treatment) %>% 
    group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
    mutate(var_comp = 1/(sigma2 + n * tau2))
  IF_ancova <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res
  est_ancova <- beta[2]
  mr_se_ancova <- sd(IF_ancova)/sqrt(m-2-p)
  mb_se_ancova <- zz$coefficients[2,2]
  
  # mix_ancova2
  mix_ancova2 <- lmer(outcome~treatment*(xdemean+Sdemean) + (1|cluster), 
                      data = sim_data %>% mutate(xdemean = x - mu_x_hat, Sdemean = S - mu_S_hat), REML = F)
  zz <- summary(mix_ancova2)
  beta <- zz$coefficients[,1]
  tau2 <- zz$varcor$cluster[1]
  sigma2 <- zz$sigma^2
  cl_res <- data.frame(res = sim_data$outcome - model.matrix(mix_ancova2) %*% beta, 
                       cluster = sim_data$cluster, 
                       a = sim_data$treatment,
                       x = sim_data$x,
                       S = sim_data$S) %>% 
    group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a), sx = sum(x), sS = sum(S)) %>%
    mutate(var_comp = 1/(sigma2 + n * tau2))
  IF_ancova2 <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res + 
    beta[5]/mean(cl_res$n) * (cl_res$sx - cl_res$n * mu_x_hat) +
    beta[6]/mean(cl_res$n) * (cl_res$sS - cl_res$n * mu_S_hat)
  est_ancova2 <- beta[2]
  mr_se_ancova2 <- sd(IF_ancova2)/sqrt(m-2-2*p)
  mb_se_ancova2 <- zz$coefficients[2,2]
  
  # individual-ancova
  indi_ancova <- lm(outcome~treatment+x+S, data = sim_data)
  zz <- summary(indi_ancova)
  beta <- zz$coefficients[,1]
  cl_res <- data.frame(res = sim_data$outcome - model.matrix(indi_ancova) %*% beta, 
                       cluster = sim_data$cluster, 
                       a = sim_data$treatment) %>% 
    group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) 
  IF_indi_ancova <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n) * cl_res$res
  est_indi_ancova <- beta[2]
  mr_se_indi_ancova <- sd(IF_indi_ancova)/sqrt(m-2-p)
  mb_se_indi_ancova <- zz$coefficients[2,2]
  
  # individual-ancova2
  indi_ancova2 <- lm(outcome~treatment*(xdemean+Sdemean),
                     data = sim_data %>% mutate(xdemean = x - mu_x_hat, Sdemean = S - mu_S_hat))
  zz <- summary(indi_ancova2)
  beta <- zz$coefficients[,1]
  cl_res <- data.frame(res = sim_data$outcome - model.matrix(indi_ancova2) %*% beta, 
                       cluster = sim_data$cluster, 
                       a = sim_data$treatment,
                       x = sim_data$x,
                       S = sim_data$S) %>% 
    group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a), sx = sum(x), sS = sum(S)) 
  IF_indi_ancova2 <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n) * cl_res$res +
    beta[5]/mean(cl_res$n) * (cl_res$sx - cl_res$n * mu_x_hat) +
    beta[6]/mean(cl_res$n) * (cl_res$sS - cl_res$n * mu_S_hat)
  est_indi_ancova2 <- beta[2]
  mr_se_indi_ancova2 <- sd(IF_indi_ancova2)/sqrt(m-2-2*p)
  mb_se_indi_ancova2 <- zz$coefficients[2,2]
  
  # cl-ancova
  cl_ancova <- lm(outcome~trt+x+S, data = cluster_level_data)
  IF_cl_ancova <- (cluster_level_data$trt - pi)/pi/(1-pi) * cl_ancova$residuals
  est_cl_ancova <- coef(cl_ancova)[2]
  mr_se_cl_ancova <- sd(IF_cl_ancova)/sqrt(m-2-p)
  mb_se_cl_ancova <- summary(cl_ancova)$coefficients[2,2]
  
  # cl-ancova2
  cl_ancova2 <- lm(outcome~trt*(xdemean+Sdemean), 
                   data = cluster_level_data %>% mutate(xdemean = x - mean(cluster_level_data$x), Sdemean = S - mean(cluster_level_data$S)))
  beta <- coef(cl_ancova2)
  IF_cl_ancova2 <- (cluster_level_data$trt - pi)/pi/(1-pi) * cl_ancova2$residuals +  
    beta["trt:xdemean"] * (cluster_level_data$x - mean(cluster_level_data$x)) +
    beta["trt:Sdemean"] * (cluster_level_data$S - mean(cluster_level_data$S))
  est_cl_ancova2 <- beta[2]
  mr_se_cl_ancova2 <- sd(IF_cl_ancova2)/sqrt(m-2-2*p)
  mb_se_cl_ancova2 <- summary(cl_ancova2)$coefficients[2,2]
  
  c(est_ancova, est_indi_ancova, est_cl_ancova, est_ancova2, est_indi_ancova2, est_cl_ancova2,
    mr_se_ancova, mr_se_indi_ancova, mr_se_cl_ancova, mr_se_ancova2, mr_se_indi_ancova2, mr_se_cl_ancova2,
    mb_se_ancova, mb_se_indi_ancova, mb_se_cl_ancova, mb_se_ancova2, mb_se_indi_ancova2, mb_se_cl_ancova2)
  
  }, error = function(e){return(rep(NA,18))})
}
stopCluster(cl)
tictoc::toc()

print("Scenario 2")
summary_result2 <- data.frame(bias = apply(sim_result[1:6,], 1, mean, na.rm = T),
                             ese = apply(sim_result[1:6,], 1, sd, na.rm = T),
                             mr_ase = apply(sim_result[7:12,], 1, mean, na.rm = T),
                             mb_ase = apply(sim_result[13:18,], 1, mean, na.rm = T), 
                             cp = map_dbl(1:6, function(j){
                               mean(abs(sim_result[j,]) - qnorm(0.975, mean=0, sd=sim_result[j+6,]) <=0, na.rm = T)
                             })) %>% mutate(re = ese[1]^2/ese^2)

rownames(summary_result2) <- c("mix-ANCOVA", "indi-ANCOVA", "cl-ANCOVA", "mix-ANCOVA2", "indi-ANCOVA2", "cl-ANCOVA2")
round(summary_result2,3)
xtable::xtable(summary_result2, digits = 2)



# Monte Carlo SE
est <- sim_result[1:6,]
se <- sim_result[7:12,]
MCSE <- cbind(apply(est, 1, sd, na.rm = T)/sqrt(n_sim),
              apply(est, 1, sd, na.rm = T)/sqrt(2 * n_sim - 1), 
              apply(se, 1, sd, na.rm = T)/sqrt(4 * n_sim * apply(se, 1, mean, na.rm = T)), 
              sqrt(summary_result2$cp * (1-summary_result2$cp)/n_sim),
              2 * (apply(se, 1, mean, na.rm = T)/mean(se[1,], na.rm = T))^2 * sqrt(1-cor(t(est), use = "complete.obs")[,1]^2)/sqrt(n_sim-1)
)
colnames(MCSE) <- c("bias", "SE", "ASE", "CP", "RE")
rownames(MCSE) <- c("mix-ANCOVA", "indi-ANCOVA", "cl-ANCOVA", "mix-ANCOVA2", "indi-ANCOVA2", "cl-ANCOVA2")
MCSE
