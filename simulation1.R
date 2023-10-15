rm(list = ls())
set.seed(123)
library(tidyverse)
library(lme4)
library(foreach)
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)


# sim-1-1 --------------------------
n_sim <- 10000
m <- 200
n_i_mean <- 8
n <- n_i_mean + 4
p <- 1
scenario <- "3" # 1, 2, or 3
# scenario 1: variance: uandj < ancova < ancova2
# scenario 2: variance: ancova < ancova2 < unadj
# sceniaro 3: variance: ancova2 < unadj < ancova

tictoc::tic()
sim_result <- foreach(k = 1:n_sim, .combine = cbind, .packages = c("tidyverse", "lme4")) %dopar% {
  sim_data <- map_dfr(1:m, function(i){
    x <- rnorm(n, sd = 2)
    if(scenario == "1"){
      pi <- 0.5
      a <- rbinom(1,size = 1, prob=pi)
      y <-  4 * a * (x- mean(x)) + rnorm(n, sd = 5)  + rnorm(1, sd = 1)
    }else if(scenario == "2"){
      pi <- 0.4
      a <- rbinom(1,size = 1, prob=pi)
      y <-  2 * x + 4 * a * (x- mean(x)) + rnorm(n, sd = 5)  + rnorm(1, sd = 1)
    } else if (scenario == "3"){
      pi <- 0.6
      a <- rbinom(1,size = 1, prob=pi)
      y <-  4 * a * (2 * x- mean(x)) + rnorm(n, sd = 5)  + rnorm(1, sd = 1)
    }
    n_i <- round(runif(1, min = n_i_mean-4, max = n_i_mean+4))
    index <- sample(1:n, size = n_i, replace = F)
    x <- x[index]
    y <- y[index]
    data.frame(outcome = y,
               cluster = as.factor(i),
               x = x,
               treatment= a)
  })
  if(scenario == "1"){ pi <- 0.5  }
  else if(scenario == "2"){pi <- 0.4  } 
  else if (scenario == "3"){pi <- 0.6   }
  
  tryCatch({
  
  # mix_unadj
  mix_unadj <- lmer(outcome~treatment + (1|cluster), data = sim_data, REML = F)
  zz <- summary(mix_unadj)
  beta <- zz$coefficients[,1]
  tau2 <- zz$varcor$cluster[1]
  sigma2 <- zz$sigma^2
  cl_res <- data.frame(res = sim_data$outcome - model.matrix(mix_unadj) %*% beta, cluster = sim_data$cluster, a = sim_data$treatment) %>% 
    group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
    mutate(var_comp = 1/(sigma2 + n * tau2))
  IF_unadj <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res
  est_unadj <- beta[2]
  mr_se_unadj <- sd(IF_unadj)/sqrt(m-2)
  # mb_se_unadj <- summary(lm(outcome~treatment, data = sim_data))$coefficients[2,2]
  # zz <- summary(mix_unadj)
  mb_se_unadj <- zz$coefficients[2,2]
  
  # mix_ancova
  mix_ancova <- lmer(outcome~treatment+x + (1|cluster), data = sim_data, REML = F)
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
  mu_x_hat <- mean(sim_data$x)
  mix_ancova2 <- lmer(outcome~treatment*xdemean + (1|cluster), data = sim_data %>% mutate(xdemean = x - mu_x_hat), REML = F)
  zz <- summary(mix_ancova2)
  beta <- zz$coefficients[,1]
  tau2 <- zz$varcor$cluster[1]
  sigma2 <- zz$sigma^2
  cl_res <- data.frame(res = sim_data$outcome - model.matrix(mix_ancova2) %*% beta, 
                       cluster = sim_data$cluster, 
                       a = sim_data$treatment,
                       x = sim_data$x) %>% 
    group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a), sx = sum(x)) %>%
    mutate(var_comp = 1/(sigma2 + n * tau2))
  IF_ancova2 <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res + 
    beta[4]/mean(cl_res$n) * (cl_res$sx - cl_res$n * mean(sim_data$x))
  est_ancova2 <- beta[2]
  mr_se_ancova2 <- sd(IF_ancova2)/sqrt(m-2-p)
  mb_se_ancova2 <- zz$coefficients[2,2]

  c(est_unadj, est_ancova, est_ancova2, 
    mr_se_unadj, mr_se_ancova, mr_se_ancova2, 
    mb_se_unadj, mb_se_ancova, mb_se_ancova2)
  
  }, error = function(e){return(rep(NA,9))})
}
stopCluster(cl)
tictoc::toc()

print("simulation-1")
summary_result1 <- data.frame(bias = apply(sim_result[1:3,], 1, mean, na.rm=T),
                             ese = apply(sim_result[1:3,], 1, sd, na.rm=T),
                             mr_ase = apply(sim_result[4:6,], 1, mean, na.rm=T),
                             mb_ase = apply(sim_result[7:9,], 1, mean, na.rm=T),
                             mr_cp = map_dbl(1:3, function(j){
                               mean(abs(sim_result[j,]) - qnorm(0.975, mean=0, sd=sim_result[j+3,]) <=0, na.rm=T)
                             }),
                             mb_cp = map_dbl(1:3, function(j){
                               mean(abs(sim_result[j,]) - qnorm(0.975, mean=0, sd=sim_result[j+6,]) <=0, na.rm=T)
                             })) %>% mutate(re = ese[1]^2/ese^2)
rownames(summary_result1) <- c("mix-unadj", "mix-ANCOVA", "mix-ANCOVA2")

xtable::xtable(summary_result1, digits = 2)

# Monte Carlo SE
est <- sim_result[1:3,]
se <- sim_result[4:6,]
MCSE <- cbind(apply(est, 1, sd, na.rm=T)/sqrt(n_sim),
              apply(est, 1, sd, na.rm=T)/sqrt(2 * n_sim - 1), 
              apply(se, 1, sd, na.rm=T)/sqrt(4 * n_sim * apply(se, 1, mean, na.rm=T)), 
              apply(sim_result[7:9,], 1, sd, na.rm=T)/sqrt(4 * n_sim * apply(sim_result[7:9,], 1, mean, na.rm=T)), 
              sqrt(summary_result1$mr_cp * (1-summary_result1$mr_cp)/n_sim),
              sqrt(summary_result1$mb_cp * (1-summary_result1$mb_cp)/n_sim),
              2 * (apply(se, 1, mean, na.rm=T)/mean(se[1,], na.rm=T))^2 * sqrt(1-cor(t(est))[,1]^2)/sqrt(n_sim-1)
)
colnames(MCSE) <- c("bias", "SE", "AmrSE", "AmbSE", "mrCP", "mbCP", "RE")
rownames(MCSE) <- c("mix-unadj", "mix-ANCOVA", "mix-ANCOVA2")
apply(MCSE, 2, max)




