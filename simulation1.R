rm(list = ls())
set.seed(123)
library(tidyverse)
library(lme4)
library(foreach)
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

n_sim <- 10000
m <- 200
n_i_mean <- 8
n <- n_i_mean + 4

tictoc::tic()
sim_result <- foreach(k = 1:n_sim, .combine = cbind, .packages = c("tidyverse", "lme4")) %dopar% {
  sim_data <- map_dfr(1:m, function(i){
    x <- rnorm(n, sd = 2)
    y <- x - mean(x) + rnorm(n, sd = 5)  + rnorm(1, sd = 1)
    n_i <- round(runif(1, min = n_i_mean-4, max = n_i_mean+4))
    index <- sample(1:n, size = n_i, replace = F)
    x <- x[index]
    y <- y[index]
    data.frame(outcome = y,
               cluster = as.factor(i),
               x = x,
               treatment= rbinom(1,size = 1, prob=0.5))
  })
  mix_unadj <- lmer(outcome~treatment + (1|cluster), data = sim_data, REML = F)
  mix_unadj_reml <- lmer(outcome~treatment + (1|cluster), data = sim_data, REML = T)
  mix_fit <- lmer(outcome~treatment+x + (1|cluster), data = sim_data, REML = F)
  mix_fit_reml <- lmer(outcome~treatment+x + (1|cluster), data = sim_data, REML = T)
  
  
  cluster_level_data <- sim_data %>% group_by(cluster) %>% 
    summarise(outcome = mean(outcome),
              x = mean(x),
              treatment = mean(treatment))
  cl_fit <- lm(outcome~treatment+x, data = cluster_level_data)
  unadjusted <- lm(outcome~treatment, data = cluster_level_data)
  
  c(summary(unadjusted)$coefficients[2,1:2],
    summary(mix_unadj)$coefficients[2,1:2],
    summary(mix_unadj_reml)$coefficients[2,1:2],
    summary(mix_fit)$coefficients[2,1:2],
    summary(mix_fit_reml)$coefficients[2,1:2],
    summary(cl_fit)$coefficients[2,1:2])
}
tictoc::toc()

print("Scenario 1")
summary_result <- data.frame(bias = apply(sim_result[c(1,3,5,7,9,11),], 1, mean),
                             se = apply(sim_result[c(1,3,5,7,9,11),], 1, sd),
                             ase = apply(sim_result[c(2,4,6,8,10,12),], 1, mean),
                             cp = map_dbl(c(1,3,5,7,9,11), function(j){
                               mean(abs(sim_result[j,]) - qnorm(0.975, mean=0, sd=sim_result[j+1,]) <=0)
                             }))
summary_result$ase[2] <- summary_result$ase[2] * sqrt(m/(m-2))
summary_result$ase[4] <- summary_result$ase[4] * sqrt(m/(m-3))
rownames(summary_result) <- c("unadj", "mix-unadj", "mix-unadj-reml", "mix-ANCOVA", "mix-ANCOVA-reml", "cl-ANCOVA")

xtable::xtable(summary_result, digits = 3)

stopCluster(cl)
