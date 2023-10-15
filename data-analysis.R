rm(list = ls())
set.seed(123)
library(tidyverse)
library(lme4)
library(readxl)
library(readstata13)

# TASSH -----------
d_outcome <- read_xlsx("Ogedegbe_2018/Observed Outcome Data.xlsx") %>% 
  dplyr::select(`Participant ID`, `SBP 12 Months`)

d_bl <- read_xlsx("Ogedegbe_2018/Clinical Data.xlsx") %>% 
  dplyr::select(`Participant ID`, Site = `Site Number`, SBP, DBP)

d_bl2 <- read_xlsx("Ogedegbe_2018/Demographic Data.xlsx") %>% 
  dplyr::select(`Participant ID`, Age, Gender)

d_trt <- read_xlsx("Ogedegbe_2018/Site Characteristic Data.xlsx") %>%
  dplyr::select(Site , treatment = `Tx or Ctrl`, Rural = `Rural/Urban`)

d <- left_join(d_trt, d_bl, by = "Site") %>%
  left_join(d_bl2, by = "Participant ID") %>%
  left_join(d_outcome, by = "Participant ID") %>%
  filter(!is.na(`SBP 12 Months`)) %>%
  mutate(Site = as.factor(Site))
d$SBP[is.na(d$SBP)] <- mean(d$SBP, na.rm = T)
d$DBP[is.na(d$DBP)] <- mean(d$DBP, na.rm = T)
d$Age[is.na(d$Age)] <- mean(d$Age, na.rm = T)
dd <- d %>% group_by(Site) %>% summarise(cluster_size = n())
d <- left_join(d, dd, by = "Site")
summary(d)
pi <- 0.5
p <- 4
m <- 32

# mix-unadj
mix_unadj <- lmer(`SBP 12 Months`~treatment + (1|Site), data = d, REML = F)
zz <- summary(mix_unadj)
beta <- zz$coefficients[,1]
tau2 <- zz$varcor$Site[1]
sigma2 <- zz$sigma^2
cl_res <- data.frame(res = d$`SBP 12 Months` - model.matrix(mix_unadj) %*% beta, 
                     cluster = d$Site, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
  mutate(var_comp = 1/(sigma2 + n * tau2))
IF_unadj <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res
est_unadj <- beta[2]
mr_se_unadj <- sd(IF_unadj)/sqrt(m-2)
mb_se_unadj <- zz$coefficients[2,2]

# mix_ancova
mix_ancova <- lmer(`SBP 12 Months`~treatment+ Rural + SBP + DBP + Age + (1|Site), data = d, REML = F)
zz <- summary(mix_ancova)
beta <- zz$coefficients[,1]
tau2 <- zz$varcor$Site[1]
sigma2 <- zz$sigma^2
cl_res <- data.frame(res = d$`SBP 12 Months` - model.matrix(mix_ancova) %*% beta, 
                     cluster = d$Site, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
  mutate(var_comp = 1/(sigma2 + n * tau2))
IF_ancova <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res
est_ancova <- beta[2]
mr_se_ancova <- sd(IF_ancova)/sqrt(m-2-p)
mb_se_ancova <- zz$coefficients[2,2]

# mix_ancova2
dc <- d %>% mutate(Ruralc = Rural - mean(Rural), SBPc = SBP - mean(SBP),
                   DBPc = DBP - mean(DBP), Agec = Age - mean(Age))
mix_ancova2 <- lmer(`SBP 12 Months`~treatment * (Ruralc + SBPc + DBPc + Agec) + (1|Site), data = dc, REML = F)
zz <- summary(mix_ancova2)
beta <- zz$coefficients[,1]
tau2 <- zz$varcor$Site[1]
sigma2 <- zz$sigma^2
cl_res <- data.frame(res = d$`SBP 12 Months` - model.matrix(mix_ancova2) %*% beta, 
                     cluster = d$Site, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
  mutate(var_comp = 1/(sigma2 + n * tau2))
beta_AX <- beta[7:10]
sX <- d[,c("Site", "Rural", "SBP", "DBP", "Age")] %>% group_by(Site) %>% summarise_all(sum) %>% .[,-1] %>% as.matrix
barX <- colMeans(d[,c("Rural", "SBP", "DBP", "Age")])
IF_ancova2 <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res + 
  1/mean(cl_res$n) * (sX %*% beta_AX - cl_res$n * c(barX %*% beta_AX))
est_ancova2 <- beta[2]
mr_se_ancova2 <- sd(IF_ancova2)/sqrt(m-2-2*p)
mb_se_ancova2 <- zz$coefficients[2,2]

# individual-ancova
indi_ancova <- lm(`SBP 12 Months`~treatment+ Rural + SBP + DBP + Age, data = d)
zz <- summary(indi_ancova)
beta <- zz$coefficients[,1]
cl_res <- data.frame(res = d$`SBP 12 Months` - model.matrix(indi_ancova) %*% beta, 
                     cluster = d$Site, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) 
IF_indi_ancova <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n) * cl_res$res
est_indi_ancova <- beta[2]
mr_se_indi_ancova <- sd(IF_indi_ancova)/sqrt(m-2-p)
mb_se_indi_ancova <- zz$coefficients[2,2]

# individual-ancova2
dc <- d %>% mutate(Ruralc = Rural - mean(Rural), SBPc = SBP - mean(SBP),
                   DBPc = DBP - mean(DBP), Agec = Age - mean(Age))
indi_ancova2 <- lm(`SBP 12 Months`~treatment * (Ruralc + SBPc + DBPc + Agec), data = dc)
zz <- summary(indi_ancova2)
beta <- zz$coefficients[,1]
cl_res <- data.frame(res = d$`SBP 12 Months` - model.matrix(indi_ancova2) %*% beta, 
                     cluster = d$Site, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) 
beta_AX <- beta[7:10]
sX <- d[,c("Site", "Rural", "SBP", "DBP", "Age")] %>% group_by(Site) %>% summarise_all(sum) %>% .[,-1] %>% as.matrix
barX <- colMeans(d[,c("Rural", "SBP", "DBP", "Age")])
IF_indi_ancova2 <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n) * cl_res$res +
  1/mean(cl_res$n) * (sX %*% beta_AX - cl_res$n * c(barX %*% beta_AX))
est_indi_ancova2 <- beta[2]
mr_se_indi_ancova2 <- sd(IF_indi_ancova2)/sqrt(m-2-2*p)
mb_se_indi_ancova2 <- zz$coefficients[2,2]

# cl-ancova
cluster_level_data <- d %>% dplyr::select(-`Participant ID`) %>% group_by(Site) %>% summarise_all(mean)
cl_ancova <- lm(`SBP 12 Months`~treatment+ Rural + SBP + DBP + Age, data = cluster_level_data)
IF_cl_ancova <- (cluster_level_data$treatment - pi)/pi/(1-pi) * cl_ancova$residuals
est_cl_ancova <- coef(cl_ancova)[2]
mr_se_cl_ancova <- sd(IF_cl_ancova)/sqrt(m-2-p)
mb_se_cl_ancova <- summary(cl_ancova)$coefficients[2,2]

# cl-ancova2
cl_data <- d %>% dplyr::select(-`Participant ID`) %>% group_by(Site) %>% summarise_all(mean)
cl_datac <- cl_data %>% mutate(Ruralc = Rural - mean(Rural), SBPc = SBP - mean(SBP),
                         DBPc = DBP - mean(DBP), Agec = Age - mean(Age))
cl_ancova2 <- lm(`SBP 12 Months`~treatment * (Ruralc + SBPc + DBPc + Agec), data = cl_datac)
beta <- coef(cl_ancova2)
betaAX <- beta[7:10]
IF_cl_ancova2 <- (cluster_level_data$treatment - pi)/pi/(1-pi) * cl_ancova2$residuals +  
  as.matrix(cl_datac[,c("Ruralc", "SBPc", "DBPc", "Agec")]) %*% beta_AX
est_cl_ancova2 <- beta[2]
mr_se_cl_ancova2 <- sd(IF_cl_ancova2)/sqrt(m-2-2*p)
mb_se_cl_ancova2 <- summary(cl_ancova2)$coefficients[2,2]

summary_table <- data.frame(est = c(est_unadj, est_ancova, est_indi_ancova, est_cl_ancova, est_ancova2, est_indi_ancova2, est_cl_ancova2),
                            sd = c(mr_se_unadj, mr_se_ancova, mr_se_indi_ancova, mr_se_cl_ancova, mr_se_ancova2, mr_se_indi_ancova2, mr_se_cl_ancova2)) %>%
  mutate(ci.lower = qnorm(0.025, est, sd),
                  ci.upper = qnorm(0.975, est, sd),
                  pvr = 1-sd^2/sd[1]^2)
xtable::xtable(summary_table, digits = 2)

# mix_unadj <- lmer(`SBP 12 Months`~treatment+ (1|Site), data = d, REML = F)
# mix_unadj_reml <- lmer(`SBP 12 Months`~treatment+ (1|Site), data = d, REML = T)
# mix_fit <- lmer(`SBP 12 Months`~treatment+ Rural + SBP + DBP + Age + (1|Site), data = d, REML = F)
# mix_reml <- lmer(`SBP 12 Months`~treatment+ Rural + SBP + DBP + Age + (1|Site), data = d, REML = T)
# cluster_level_data <- d %>% dplyr::select(-`Participant ID`) %>% group_by(Site) %>% summarise_all(mean)
# unadjusted <- lm(`SBP 12 Months`~treatment, data = cluster_level_data)
# cl_fit <- lm(`SBP 12 Months`~treatment+ Rural + SBP + DBP + Age, data = cluster_level_data)

# summary_table <- data.frame(est = c(summary(mix_unadj)$coefficients[2,1], 
#                                     summary(unadjusted)$coefficients[2,1], 
#                                     summary(mix_unadj_reml)$coefficients[2,1], 
#                                     summary(mix_fit)$coefficients[2,1], 
#                                     summary(mix_reml)$coefficients[2,1], 
#                                     summary(cl_fit)$coefficients[2,1]),
#                             sd = c(summary(mix_unadj)$coefficients[2,2] * sqrt(32/30), 
#                                    summary(unadjusted)$coefficients[2,2], 
#                                    summary(mix_unadj_reml)$coefficients[2,2], 
#                                    summary(mix_fit)$coefficients[2,2] * sqrt(32/26), 
#                                    summary(mix_reml)$coefficients[2,2],
#                                    summary(cl_fit)$coefficients[2,2])) %>%
#   mutate(ci.lower = qnorm(0.025, est, sd),
#          ci.upper = qnorm(0.975, est, sd),
#          pvr = 1-sd^2/sd[1]^2)
# rownames(summary_table) <- c("mix-unadj", "unadj", "mix-unadj-reml", "mix-ANCOVA", "mix-ANCOVA-reml", "cl-ANCOVA")

# xtable::xtable(summary_table, digits = 3)

# IECDZ ------------
d <- read.csv("Rockers_2018/Rockers.csv") %>% 
  dplyr::select(cluster, treatment, age_months, sex, haz06, haz_bl, waz_bl, whz_bl, SB_motor_bl) %>%
  filter(!is.na(haz06)) %>%
  mutate(cluster = as.factor(cluster))
d$haz_bl[is.na(d$haz_bl)] <- mean(d$haz_bl, na.rm = T)
d$waz_bl[is.na(d$waz_bl)] <- mean(d$waz_bl, na.rm = T)
d$whz_bl[is.na(d$whz_bl)] <- mean(d$whz_bl, na.rm = T)
d$SB_motor_bl[is.na(d$SB_motor_bl )] <- mean(d$SB_motor_bl, na.rm = T)
dd <- d %>% group_by(cluster) %>% summarise(cluster_size = n())
d <- left_join(d, dd, by = "cluster")
summary(d)
pi <- 0.5
p <- 3
m <- 30

# mix-unadj
mix_unadj <- lmer(haz06~treatment + (1|cluster), data = d, REML = F)
zz <- summary(mix_unadj)
beta <- zz$coefficients[,1]
tau2 <- zz$varcor$cluster[1]
sigma2 <- zz$sigma^2
cl_res <- data.frame(res = d$haz06 - model.matrix(mix_unadj) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
  mutate(var_comp = 1/(sigma2 + n * tau2))
IF_unadj <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res
est_unadj <- beta[2]
mr_se_unadj <- sd(IF_unadj)/sqrt(m-2)
mb_se_unadj <- zz$coefficients[2,2]

# mix_ancova
mix_ancova <- lmer(haz06~treatment+ haz_bl+SB_motor_bl+age_months + (1|cluster), data = d, REML = F)
zz <- summary(mix_ancova)
beta <- zz$coefficients[,1]
tau2 <- zz$varcor$cluster[1]
sigma2 <- zz$sigma^2
cl_res <- data.frame(res = d$haz06 - model.matrix(mix_ancova) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
  mutate(var_comp = 1/(sigma2 + n * tau2))
IF_ancova <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res
est_ancova <- beta[2]
mr_se_ancova <- sd(IF_ancova)/sqrt(m-2-p)
mb_se_ancova <- zz$coefficients[2,2]

# mix_ancova2
dc <- d %>% mutate(haz_blc = haz_bl - mean(haz_bl), SB_motor_blc = SB_motor_bl - mean(SB_motor_bl),
                   age_monthsc = age_months - mean(age_months))
mix_ancova2 <- lmer(haz06~treatment * (haz_blc+SB_motor_blc+age_monthsc) + (1|cluster), data = dc, REML = F)
zz <- summary(mix_ancova2)
beta <- zz$coefficients[,1]
tau2 <- zz$varcor$cluster[1]
sigma2 <- zz$sigma^2
cl_res <- data.frame(res = d$haz06 - model.matrix(mix_ancova2) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
  mutate(var_comp = 1/(sigma2 + n * tau2))
beta_AX <- beta[6:8]
sX <- d[,c("cluster", "haz_bl", "SB_motor_bl", "age_months")] %>% group_by(cluster) %>% summarise_all(sum) %>% .[,-1] %>% as.matrix
barX <- colMeans(d[,c("haz_bl", "SB_motor_bl", "age_months")])
IF_ancova2 <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res + 
  1/mean(cl_res$n) * (sX %*% beta_AX - cl_res$n * c(barX %*% beta_AX))
est_ancova2 <- beta[2]
mr_se_ancova2 <- sd(IF_ancova2)/sqrt(m-2-2*p)
mb_se_ancova2 <- zz$coefficients[2,2]

# individual-ancova
indi_ancova <- lm(haz06~treatment+ haz_bl+SB_motor_bl+age_months, data = d)
zz <- summary(indi_ancova)
beta <- zz$coefficients[,1]
cl_res <- data.frame(res = d$haz06 - model.matrix(indi_ancova) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) 
IF_indi_ancova <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n) * cl_res$res
est_indi_ancova <- beta[2]
mr_se_indi_ancova <- sd(IF_indi_ancova)/sqrt(m-2-p)
mb_se_indi_ancova <- zz$coefficients[2,2]

# individual-ancova2
dc <- d %>% mutate(haz_blc = haz_bl - mean(haz_bl), SB_motor_blc = SB_motor_bl - mean(SB_motor_bl),
                   age_monthsc = age_months - mean(age_months))
indi_ancova2 <- lm(haz06~treatment * (haz_blc+SB_motor_blc+age_monthsc), data = dc)
zz <- summary(indi_ancova2)
beta <- zz$coefficients[,1]
cl_res <- data.frame(res = d$haz06 - model.matrix(indi_ancova2) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) 
beta_AX <- beta[6:8]
sX <- d[,c("cluster", "haz_bl", "SB_motor_bl", "age_months")] %>% group_by(cluster) %>% summarise_all(sum) %>% .[,-1] %>% as.matrix
barX <- colMeans(d[,c("haz_bl", "SB_motor_bl", "age_months")])
IF_indi_ancova2 <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n) * cl_res$res +
  1/mean(cl_res$n) * (sX %*% beta_AX - cl_res$n * c(barX %*% beta_AX))
est_indi_ancova2 <- beta[2]
mr_se_indi_ancova2 <- sd(IF_indi_ancova2)/sqrt(m-2-2*p)
mb_se_indi_ancova2 <- zz$coefficients[2,2]

# cl-ancova
cluster_level_data <- d %>% group_by(cluster) %>% summarise_all(mean)
cl_ancova <- lm(haz06~treatment+ haz_bl+SB_motor_bl+age_months, data = cluster_level_data)
IF_cl_ancova <- (cluster_level_data$treatment - pi)/pi/(1-pi) * cl_ancova$residuals
est_cl_ancova <- coef(cl_ancova)[2]
mr_se_cl_ancova <- sd(IF_cl_ancova)/sqrt(m-2-p)
mb_se_cl_ancova <- summary(cl_ancova)$coefficients[2,2]

# cl-ancova2
cl_data <- d %>% group_by(cluster) %>% summarise_all(mean)
cl_datac <- cl_data %>% mutate(haz_blc = haz_bl - mean(haz_bl), SB_motor_blc = SB_motor_bl - mean(SB_motor_bl),
                               age_monthsc = age_months - mean(age_months))
cl_ancova2 <- lm(haz06~treatment * (haz_blc+SB_motor_blc+age_monthsc), data = cl_datac)
beta <- coef(cl_ancova2)
beta_AX <- beta[6:8]
IF_cl_ancova2 <- (cl_data$treatment - pi)/pi/(1-pi) * cl_ancova2$residuals +  
  as.matrix(cl_datac[,c("haz_bl", "SB_motor_bl", "age_months")]) %*% beta_AX
est_cl_ancova2 <- beta[2]
mr_se_cl_ancova2 <- sd(IF_cl_ancova2)/sqrt(m-2-2*p)
mb_se_cl_ancova2 <- summary(cl_ancova2)$coefficients[2,2]

summary_table <- data.frame(est = c(est_unadj, est_ancova, est_indi_ancova, est_cl_ancova, est_ancova2, est_indi_ancova2, est_cl_ancova2),
                            sd = c(mr_se_unadj, mr_se_ancova, mr_se_indi_ancova, mr_se_cl_ancova, mr_se_ancova2, mr_se_indi_ancova2, mr_se_cl_ancova2)) %>%
  mutate(ci.lower = qnorm(0.025, est, sd),
         ci.upper = qnorm(0.975, est, sd),
         pvr = 1-sd^2/sd[1]^2)
rownames(summary_table) <- c("mix-unadj", "mix-ANCOVA", "indi-ancova", "cl-ancova", "mix-ANCOVA2","indi-ancova2", "cl-ancova2")
xtable::xtable(summary_table, digits = 2)


# mix_unadj <- lmer(haz06~treatment + (1|cluster), data = d, REML = F)
# mix_unadj_reml <- lmer(haz06~treatment + (1|cluster), data = d, REML = T)
# mix_fit <- lmer(haz06~treatment+haz_bl+SB_motor_bl+age_months + (1|cluster), data = d, REML = F)
# mix_reml <- lmer(haz06~treatment+haz_bl+SB_motor_bl+age_months + (1|cluster), data = d, REML = T)
# cluster_level_data <- d %>% group_by(cluster) %>% summarise_all(mean)
# unadjusted <- lm(haz06~treatment, data = cluster_level_data)
# cl_fit <- lm(haz06~treatment+haz_bl+SB_motor_bl+age_months, data = cluster_level_data)
# 
# 
# summary_table <- data.frame(est = c(summary(mix_unadj)$coefficients[2,1], 
#                                     summary(unadjusted)$coefficients[2,1], 
#                                     summary(mix_unadj_reml)$coefficients[2,1], 
#                                     summary(mix_fit)$coefficients[2,1], 
#                                     summary(mix_reml)$coefficients[2,1], 
#                                     summary(cl_fit)$coefficients[2,1]),
#                             sd = c(summary(mix_unadj)$coefficients[2,2] * sqrt(30/28), 
#                                    summary(unadjusted)$coefficients[2,2], 
#                                    summary(mix_unadj_reml)$coefficients[2,2], 
#                                    summary(mix_fit)$coefficients[2,2] * sqrt(30/25), 
#                                    summary(mix_reml)$coefficients[2,2],
#                                    summary(cl_fit)$coefficients[2,2])) %>%
#   mutate(ci.lower = qnorm(0.025, est, sd),
#          ci.upper = qnorm(0.975, est, sd),
#          pvr = 1-sd^2/sd[1]^2)
# rownames(summary_table) <- c("mix-unadj", "unadj", "mix-unadj-reml", "mix-ANCOVA", "mix-ANCOVA-reml", "cl-ANCOVA")
# 
# xtable::xtable(summary_table, digits = 3)

# work-family ----------
d <- read.csv("36158-0001-Data.csv") %>%
  dplyr::select(ADMINLINK, WAVE, EMPLOYEE, STUDYGROUP, CONDITION, RMZFN, RMZEMP, SCWM_CWH) %>%
  mutate(SCWM_CWH = ifelse(SCWM_CWH >=0, SCWM_CWH, NA)) %>%
  filter(EMPLOYEE == 1 & WAVE %in% c(1,2)) %>%
  mutate(WAVE = ifelse(WAVE==1, "baseline", "six_month")) %>%
  pivot_wider(names_from = WAVE, values_from = SCWM_CWH) %>%
  mutate(treatment = ifelse(CONDITION==1, 1, 0), cluster = as.factor(STUDYGROUP)) %>%
  dplyr::select(treatment, cluster, RMZFN, RMZEMP, baseline, six_month) %>%
  filter(!is.na(six_month))
d$baseline[is.na(d$baseline)] <- mean(d$baseline, na.rm = T)
d$RMZFN[is.na(d$RMZFN)] <- median(d$RMZFN, na.rm = T)
d$RMZEMP[is.na(d$RMZEMP)] <- mean(d$RMZEMP, na.rm = T)
summary(d)
pi <- 0.5
p <- 3
m <- 56

# mix-unadj
mix_unadj <- lmer(six_month~treatment + (1|cluster), data = d, REML = F)
zz <- summary(mix_unadj)
beta <- zz$coefficients[,1]
tau2 <- zz$varcor$cluster[1]
sigma2 <- zz$sigma^2
cl_res <- data.frame(res = d$six_month - model.matrix(mix_unadj) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
  mutate(var_comp = 1/(sigma2 + n * tau2))
IF_unadj <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res
est_unadj <- beta[2]
mr_se_unadj <- sd(IF_unadj)/sqrt(m-2)
mb_se_unadj <- zz$coefficients[2,2]

# mix_ancova
mix_ancova <- lmer(six_month~treatment+ RMZFN + baseline + (1|cluster), data = d, REML = F)
zz <- summary(mix_ancova)
beta <- zz$coefficients[,1]
tau2 <- zz$varcor$cluster[1]
sigma2 <- zz$sigma^2
cl_res <- data.frame(res = d$six_month - model.matrix(mix_ancova) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
  mutate(var_comp = 1/(sigma2 + n * tau2))
IF_ancova <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res
est_ancova <- beta[2]
mr_se_ancova <- sd(IF_ancova)/sqrt(m-2-p)
mb_se_ancova <- zz$coefficients[2,2]

# mix_ancova2
dc <- d %>% mutate(RMZFNc = RMZFN - mean(RMZFN), #RMZEMPc = RMZEMP - mean(RMZEMP),
                   baselinec = baseline - mean(baseline))
mix_ancova2 <- lmer(six_month~treatment * (RMZFNc+baselinec) + (1|cluster), data = dc, REML = F)
zz <- summary(mix_ancova2)
beta <- zz$coefficients[,1]
tau2 <- zz$varcor$cluster[1]
sigma2 <- zz$sigma^2
cl_res <- data.frame(res = d$six_month - model.matrix(mix_ancova2) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) %>%
  mutate(var_comp = 1/(sigma2 + n * tau2))
beta_AX <- beta[5:6]
sX <- d[,c("cluster", "RMZFN","baseline")] %>% group_by(cluster) %>% summarise_all(sum) %>% .[,-1] %>% as.matrix
barX <- colMeans(d[,c("RMZFN", "baseline")])
IF_ancova2 <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n * cl_res$var_comp) * cl_res$var_comp * cl_res$res + 
  1/mean(cl_res$n) * (sX %*% beta_AX - cl_res$n * c(barX %*% beta_AX))
est_ancova2 <- beta[2]
mr_se_ancova2 <- sd(IF_ancova2)/sqrt(m-2-2*p)
mb_se_ancova2 <- zz$coefficients[2,2]

# individual-ancova
indi_ancova <- lm(six_month~treatment+ RMZFN + baseline, data = d)
zz <- summary(indi_ancova)
beta <- zz$coefficients[,1]
cl_res <- data.frame(res = d$six_month - model.matrix(indi_ancova) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) 
IF_indi_ancova <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n) * cl_res$res
est_indi_ancova <- beta[2]
mr_se_indi_ancova <- sd(IF_indi_ancova)/sqrt(m-2-p)
mb_se_indi_ancova <- zz$coefficients[2,2]

# individual-ancova2
dc <- d %>% mutate(RMZFNc = RMZFN - mean(RMZFN), #RMZEMPc = RMZEMP - mean(RMZEMP),
                   baselinec = baseline - mean(baseline))
indi_ancova2 <- lm(six_month~treatment * (RMZFNc+baselinec), data = dc)
zz <- summary(indi_ancova2)
beta <- zz$coefficients[,1]
cl_res <- data.frame(res = d$six_month - model.matrix(indi_ancova2) %*% beta, 
                     cluster = d$cluster, 
                     a = d$treatment) %>% 
  group_by(cluster) %>% summarise(n = n(), res = sum(res), trt = mean(a)) 
beta_AX <- beta[5:6]
sX <- d[,c("cluster", "RMZFN", "baseline")] %>% group_by(cluster) %>% summarise_all(sum) %>% .[,-1] %>% as.matrix
barX <- colMeans(d[,c("RMZFN", "baseline")])
IF_indi_ancova2 <- (cl_res$trt - pi)/pi/(1-pi)/mean(cl_res$n) * cl_res$res +
  1/mean(cl_res$n) * (sX %*% beta_AX - cl_res$n * c(barX %*% beta_AX))
est_indi_ancova2 <- beta[2]
mr_se_indi_ancova2 <- sd(IF_indi_ancova2)/sqrt(m-2-2*p)
mb_se_indi_ancova2 <- zz$coefficients[2,2]

# cl-ancova
cluster_level_data <- d %>% group_by(cluster) %>% summarise_all(mean)
cl_ancova <- lm(six_month~treatment+ RMZFN + baseline, data = cluster_level_data)
IF_cl_ancova <- (cluster_level_data$treatment - pi)/pi/(1-pi) * cl_ancova$residuals
est_cl_ancova <- coef(cl_ancova)[2]
mr_se_cl_ancova <- sd(IF_cl_ancova)/sqrt(m-2-p)
mb_se_cl_ancova <- summary(cl_ancova)$coefficients[2,2]

# cl-ancova2
cl_data <- d %>% group_by(cluster) %>% summarise_all(mean)
cl_datac <- cl_data %>% mutate(RMZFNc = RMZFN - mean(RMZFN), #RMZEMPc = RMZEMP - mean(RMZEMP),
                               baselinec = baseline - mean(baseline))
cl_ancova2 <- lm(six_month~treatment * (RMZFNc+baselinec), data = cl_datac)
beta <- coef(cl_ancova2)
beta_AX <- beta[5:6]
IF_cl_ancova2 <- (cl_data$treatment - pi)/pi/(1-pi) * cl_ancova2$residuals +  
  as.matrix(cl_datac[,c("RMZFN", "baseline")]) %*% beta_AX
est_cl_ancova2 <- beta[2]
mr_se_cl_ancova2 <- sd(IF_cl_ancova2)/sqrt(m-2-2*p)
mb_se_cl_ancova2 <- summary(cl_ancova2)$coefficients[2,2]

summary_table <- data.frame(est = c(est_unadj, est_ancova, est_indi_ancova, est_cl_ancova, est_ancova2, est_indi_ancova2, est_cl_ancova2),
                            sd = c(mr_se_unadj, mr_se_ancova, mr_se_indi_ancova, mr_se_cl_ancova, mr_se_ancova2, mr_se_indi_ancova2, mr_se_cl_ancova2)) %>%
  mutate(ci.lower = qnorm(0.025, est, sd),
         ci.upper = qnorm(0.975, est, sd),
         pvr = 1-sd^2/sd[1]^2)
rownames(summary_table) <- c("mix-unadj", "mix-ANCOVA", "indi-ancova", "cl-ancova", "mix-ANCOVA2","indi-ancova2", "cl-ancova2")
xtable::xtable(summary_table, digits = 2)


# mix_unadj <- lmer(six_month ~ treatment +(1|cluster), data = d, REML = F)
# mix_unadj_reml <- lmer(six_month ~ treatment +(1|cluster), data = d, REML = T)
# mix_fit <- lmer(six_month ~ treatment + RMZFN + RMZEMP + baseline +(1|cluster), data = d, REML = F)
# mix_reml <- lmer(six_month ~ treatment + RMZFN + RMZEMP + baseline +(1|cluster), data = d, REML = T)
# cluster_level_data <- d %>% group_by(cluster) %>% summarise_all(mean)
# unadjusted <- lm(six_month ~ treatment, data = cluster_level_data)
# cl_fit <- lm(six_month ~ treatment + RMZFN + RMZEMP + baseline, data = cluster_level_data)
# summary_table <- data.frame(est = c(summary(mix_unadj)$coefficients[2,1], 
#                                     summary(unadjusted)$coefficients[2,1], 
#                                     summary(mix_unadj_reml)$coefficients[2,1], 
#                                     summary(mix_fit)$coefficients[2,1], 
#                                     summary(mix_reml)$coefficients[2,1], 
#                                     summary(cl_fit)$coefficients[2,1]),
#                             sd = c(summary(mix_unadj)$coefficients[2,2] * sqrt(56/54), 
#                                    summary(unadjusted)$coefficients[2,2], 
#                                    summary(mix_unadj_reml)$coefficients[2,2], 
#                                    summary(mix_fit)$coefficients[2,2] * sqrt(56/51), 
#                                    summary(mix_reml)$coefficients[2,2],
#                                    summary(cl_fit)$coefficients[2,2])) %>%
#   mutate(ci.lower = qnorm(0.025, est, sd),
#          ci.upper = qnorm(0.975, est, sd),
#          pvr = 1-sd^2/sd[1]^2)
# rownames(summary_table) <- c("mix-unadj", "unadj", "mix-unadj-reml", "mix-ANCOVA", "mix-ANCOVA-reml", "cl-ANCOVA")
# 
# xtable::xtable(summary_table, digits = 2)

