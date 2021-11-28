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


mix_unadj <- lmer(`SBP 12 Months`~treatment+ (1|Site), data = d, REML = F)
mix_unadj_reml <- lmer(`SBP 12 Months`~treatment+ (1|Site), data = d, REML = T)
mix_fit <- lmer(`SBP 12 Months`~treatment+ Rural + SBP + DBP + Age + (1|Site), data = d, REML = F)
mix_reml <- lmer(`SBP 12 Months`~treatment+ Rural + SBP + DBP + Age + (1|Site), data = d, REML = T)
cluster_level_data <- d %>% dplyr::select(-`Participant ID`) %>% group_by(Site) %>% summarise_all(mean)
unadjusted <- lm(`SBP 12 Months`~treatment, data = cluster_level_data)
cl_fit <- lm(`SBP 12 Months`~treatment+ Rural + SBP + DBP + Age, data = cluster_level_data)

summary_table <- data.frame(est = c(summary(mix_unadj)$coefficients[2,1], 
                                    summary(unadjusted)$coefficients[2,1], 
                                    summary(mix_unadj_reml)$coefficients[2,1], 
                                    summary(mix_fit)$coefficients[2,1], 
                                    summary(mix_reml)$coefficients[2,1], 
                                    summary(cl_fit)$coefficients[2,1]),
                            sd = c(summary(mix_unadj)$coefficients[2,2] * sqrt(32/30), 
                                   summary(unadjusted)$coefficients[2,2], 
                                   summary(mix_unadj_reml)$coefficients[2,2], 
                                   summary(mix_fit)$coefficients[2,2] * sqrt(32/26), 
                                   summary(mix_reml)$coefficients[2,2],
                                   summary(cl_fit)$coefficients[2,2])) %>%
  mutate(ci.lower = qnorm(0.025, est, sd),
         ci.upper = qnorm(0.975, est, sd),
         pvr = 1-sd^2/sd[1]^2)
rownames(summary_table) <- c("mix-unadj", "unadj", "mix-unadj-reml", "mix-ANCOVA", "mix-ANCOVA-reml", "cl-ANCOVA")

xtable::xtable(summary_table, digits = 3)

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

mix_unadj <- lmer(haz06~treatment + (1|cluster), data = d, REML = F)
mix_unadj_reml <- lmer(haz06~treatment + (1|cluster), data = d, REML = T)
mix_fit <- lmer(haz06~treatment+haz_bl+SB_motor_bl+age_months + (1|cluster), data = d, REML = F)
mix_reml <- lmer(haz06~treatment+haz_bl+SB_motor_bl+age_months + (1|cluster), data = d, REML = T)
cluster_level_data <- d %>% group_by(cluster) %>% summarise_all(mean)
unadjusted <- lm(haz06~treatment, data = cluster_level_data)
cl_fit <- lm(haz06~treatment+haz_bl+SB_motor_bl+age_months, data = cluster_level_data)


summary_table <- data.frame(est = c(summary(mix_unadj)$coefficients[2,1], 
                                    summary(unadjusted)$coefficients[2,1], 
                                    summary(mix_unadj_reml)$coefficients[2,1], 
                                    summary(mix_fit)$coefficients[2,1], 
                                    summary(mix_reml)$coefficients[2,1], 
                                    summary(cl_fit)$coefficients[2,1]),
                            sd = c(summary(mix_unadj)$coefficients[2,2] * sqrt(30/28), 
                                   summary(unadjusted)$coefficients[2,2], 
                                   summary(mix_unadj_reml)$coefficients[2,2], 
                                   summary(mix_fit)$coefficients[2,2] * sqrt(30/25), 
                                   summary(mix_reml)$coefficients[2,2],
                                   summary(cl_fit)$coefficients[2,2])) %>%
  mutate(ci.lower = qnorm(0.025, est, sd),
         ci.upper = qnorm(0.975, est, sd),
         pvr = 1-sd^2/sd[1]^2)
rownames(summary_table) <- c("mix-unadj", "unadj", "mix-unadj-reml", "mix-ANCOVA", "mix-ANCOVA-reml", "cl-ANCOVA")

xtable::xtable(summary_table, digits = 3)

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

mix_unadj <- lmer(six_month ~ treatment +(1|cluster), data = d, REML = F)
mix_unadj_reml <- lmer(six_month ~ treatment +(1|cluster), data = d, REML = T)
mix_fit <- lmer(six_month ~ treatment + RMZFN + RMZEMP + baseline +(1|cluster), data = d, REML = F)
mix_reml <- lmer(six_month ~ treatment + RMZFN + RMZEMP + baseline +(1|cluster), data = d, REML = T)
cluster_level_data <- d %>% group_by(cluster) %>% summarise_all(mean)
unadjusted <- lm(six_month ~ treatment, data = cluster_level_data)
cl_fit <- lm(six_month ~ treatment + RMZFN + RMZEMP + baseline, data = cluster_level_data)
summary_table <- data.frame(est = c(summary(mix_unadj)$coefficients[2,1], 
                                    summary(unadjusted)$coefficients[2,1], 
                                    summary(mix_unadj_reml)$coefficients[2,1], 
                                    summary(mix_fit)$coefficients[2,1], 
                                    summary(mix_reml)$coefficients[2,1], 
                                    summary(cl_fit)$coefficients[2,1]),
                            sd = c(summary(mix_unadj)$coefficients[2,2] * sqrt(56/54), 
                                   summary(unadjusted)$coefficients[2,2], 
                                   summary(mix_unadj_reml)$coefficients[2,2], 
                                   summary(mix_fit)$coefficients[2,2] * sqrt(56/51), 
                                   summary(mix_reml)$coefficients[2,2],
                                   summary(cl_fit)$coefficients[2,2])) %>%
  mutate(ci.lower = qnorm(0.025, est, sd),
         ci.upper = qnorm(0.975, est, sd),
         pvr = 1-sd^2/sd[1]^2)
rownames(summary_table) <- c("mix-unadj", "unadj", "mix-unadj-reml", "mix-ANCOVA", "mix-ANCOVA-reml", "cl-ANCOVA")

xtable::xtable(summary_table, digits = 2)

