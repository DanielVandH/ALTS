################################################################################
## Load functions and libraries
################################################################################
source("Simulation_Functions.R")
source("ALTS_Functions.R")
source("Case_Study_Functions.R")
library(lemon)
library(EnvStats)
library(L1pack)
library(caret)
library(formula.tools)
library(pracma)
library(tidyverse)
library(magrittr)
library(eeptools)
library(dplyr)
library(ggpubr)
library(MLmetrics)
library(scales)
library(ggplot2)
library(quantreg)
library(MASS)
library(viridis)
library(rlmDataDriven)

plot_aes = function(p, size = 20) {
  q = p + theme(
    plot.title = element_text(size = size),
    axis.title.x = element_text(size = size),
    axis.title.y = element_text(size = size),
    axis.text.x = element_text(size = size),
    axis.text.y = element_text(size = size),
    legend.title = element_text(size = size),
    legend.text = element_text(size = size),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    )
  ) + grids()
  theme_set(theme_gray(base_size = size))
  return(q) 
}

################################################################################
## Load all the simulation study data 
################################################################################
load("data/Simulation/SimulationStudy_Fixedq.RData") # results_fixed_q
load("data/Simulation/SimulationStudy_Varyingq.RData") # results_varying_q
load("data/Simulation/SimulationStudy_ExtraIterate.RData") # extra_iterate

p_1 = seq(0.5, 1.0, by = 0.05)
q_1 = 1.35
p_2 = seq(0.5, 1.0, by = 0.05)
q_2 = seq(1.0, 1.5, by = 0.07)
sigma_1 = 0.1
sigma_2 = 1.3

################################################################################
## Process the above results into a format appropriate for plotting
################################################################################

## results_fixed_q
new_q_fixed_q = results_fixed_q$New$q %>% t
new_p_fixed_q = results_fixed_q$New$p %>% t
new_sigma_fixed_q = results_fixed_q$New$sigma %>% t
new_rmape_fixed_q = results_fixed_q$New$rmape %>% t
bacher_p_fixed_q = results_fixed_q$Bacher$p %>% t
bacher_sigma_fixed_q = results_fixed_q$Bacher$sigma %>% t
bacher_rmape_fixed_q = results_fixed_q$Bacher$rmape %>% t
new_p_fixed_q %<>% as.data.frame %>% mutate(p = p_1) 
names(new_p_fixed_q)[1] = "phat"
new_sigma_fixed_q %<>% as.data.frame %>% mutate(p = p_1)
names(new_sigma_fixed_q)[1] = "sigmahat"
new_rmape_fixed_q %<>% as.data.frame %>% mutate(p = p_1) 
names(new_rmape_fixed_q)[1] = "mape"
new_q_fixed_q %<>% as.data.frame %>% mutate(p = p_1) 
names(new_q_fixed_q)[1] = "qhat"
bacher_p_fixed_q %<>% as.data.frame %>% mutate(p = p_1)
names(bacher_p_fixed_q)[1] = "phat"
bacher_sigma_fixed_q %<>% as.data.frame %>% mutate(p = p_1)
names(bacher_sigma_fixed_q)[1] = "sigmahat"
bacher_rmape_fixed_q %<>% as.data.frame %>% mutate(p = p_1)
names(bacher_rmape_fixed_q)[1] = "mape"

## results_varying_q
new_p_varying_q = results_varying_q$New$p %>% t 
new_sigma_varying_q = results_varying_q$New$sigma %>% t 
new_rmape_varying_q = results_varying_q$New$rmape %>% t
bacher_p_varying_q = results_varying_q$Bacher$p %>% t
bacher_sigma_varying_q = results_varying_q$Bacher$sigma %>% t
bacher_rmape_varying_q = results_varying_q$Bacher$rmape %>% t
new_p_varying_q %<>% as.data.frame %>% mutate(p = p_2) %>%
  pivot_longer(cols = 1:length(q_2), names_to = "q", values_to = "phat")
new_sigma_varying_q %<>% as.data.frame %>% mutate(p = p_2) %>%
  pivot_longer(cols = 1:length(q_2), names_to = "q", values_to = "sigmahat")
new_rmape_varying_q %<>% as.data.frame %>% mutate(p = p_2) %>%
  pivot_longer(cols = 1:length(q_2), names_to = "q", values_to = "mape")
bacher_p_varying_q %<>% as.data.frame %>% mutate(p = p_2)
names(bacher_p_varying_q)[1] = "phat"
bacher_sigma_varying_q %<>% as.data.frame %>% mutate(p = p_2)
names(bacher_sigma_varying_q)[1] = "sigmahat"
bacher_rmape_varying_q %<>% as.data.frame %>% mutate(p = p_2)
names(bacher_rmape_varying_q)[1] = "mape"

## extra_iterate

# Below is how we obtain the extra_iterate variable loaded above
  #set.seed(5536066)
  #extra_iterate = sim_fnc_1(beta, X, p, sigma_1, sigma_2, include_bacher = TRUE)
  #extra_iterate$Bacher$p = extra_iterate$Bacher$p[2:length(extra_iterate$Bacher$p)]
  #extra_iterate$Bacher$sigma = extra_iterate$Bacher$sigma[2:length(extra_iterate$Bacher$sigma)]
  #extra_iterate$New$p = extra_iterate$New$p[2:length(extra_iterate$New$p)]
  #extra_iterate$New$sigma = extra_iterate$New$sigma[2:length(extra_iterate$New$sigma)]
bacher_iterate = data.frame(Iteration = 1:length(extra_iterate$Bacher$p), p = extra_iterate$Bacher$p, sigma = extra_iterate$Bacher$sigma, Method = "Bacher")
new_iterate = data.frame(Iteration = 1:length(extra_iterate$New$p), p = extra_iterate$New$p, sigma = extra_iterate$New$sigma, Method = "New")
all_iterate = rbind(bacher_iterate, new_iterate)

################################################################################
## Figure 1: Simulation study results with fixed q.
################################################################################

## (a): MAPE values
fig_1_a = ggplot(new_rmape_fixed_q, aes(p, mape)) %>% plot_aes() +
  geom_line(color = "black", size = 2) +
  labs(y = "MAPE") +
  geom_line(data = bacher_rmape_fixed_q, mapping = aes(p, mape), size = 2, color = "red") +
  ylim(c(0.005, 0.03)) +
  scale_y_continuous(labels=scales::percent) +
  ggtitle("(a): MAPE values")

## (b): p iterates
fig_1_b = ggplot(all_iterate, aes(Iteration, p, color = Method)) %>% plot_aes() +
  geom_line(size = 2) +
  labs(y = expression(hat(p)), x = "Iteration") +
  geom_hline(yintercept = 0.75, color = "blue", linetype = "dashed") +
  scale_color_manual(values = c("red", "black")) +
  ggtitle(expression("(b): "*hat(p)*" iterates"))

## (c): sigma iterates 
fig_1_c = ggplot(all_iterate, aes(Iteration, sigma, color = Method)) %>% plot_aes() +
  geom_line(size = 2) +
  labs(y = expression(hat(sigma)), x = "Iteration") +
  geom_hline(yintercept = sigma_1, color = "blue", linetype = "dashed") +
  scale_color_manual(values = c("red", "black")) +
  ggtitle(expression("(c): "*hat(sigma)*" iterates"))

## Combine the plots 
fig_1 = ggarrange(fig_1_a, fig_1_b, fig_1_c, 
                  common.legend = TRUE,
                  nrow = 1, ncol = 3,
                  align = "h")
ggsave("figures/Figure1_SimulationStudy_Fixedq.pdf", fig_1, 
       device = "pdf", width = 13.04, height = 4.68)

################################################################################
## Figure C1: Simulation study results with varying q.
################################################################################

## (a): Estimating p
fig_c1_a = ggplot(new_p_varying_q, aes(p, phat, color = q)) %>% plot_aes() +
  geom_line(size = 2) +
  scale_color_discrete(labels = q_2) +
  labs(y = expression(hat(p))) +
  geom_abline(slope = 1, size = 2, alpha = 0.2) +
  geom_line(data = bacher_p_varying_q, aes(p, phat), color = "red", linetype = "dashed", size = 2) +
  xlim(c(0.49, 1)) +
  ylim(c(0.49, 1)) +
  ggtitle("(a): Estimating p")

## (b): Estimating sigma
fig_c1_b = ggplot(new_sigma_varying_q, aes(p, sigmahat, color = q)) %>% plot_aes() +
  geom_line(size = 2) +
  scale_color_discrete(labels = q_2) +
  labs(y = expression(hat(sigma))) +
  geom_hline(yintercept = sigma_1, size = 2) +
  geom_line(data = bacher_sigma_varying_q, aes(p, sigmahat), color = "red", linetype = "dashed", size = 2) +
  ylim(c(0.09, 0.25)) +
  ggtitle(expression("(b): Estimating "*sigma))

## (c): MAPE values 
fig_c1_c = ggplot(new_rmape_varying_q, aes(p, mape, color = q)) %>% plot_aes() +
  geom_line(size = 2) +
  scale_color_discrete(labels = q_2) +
  labs(y = "MAPE") + 
  geom_line(data = bacher_rmape_varying_q, aes(p, mape), color = "red", linetype = "dashed", size = 2) +
  ylim(c(0, 0.025)) +
  scale_y_continuous(labels=scales::percent) + 
  ggtitle("(c): MAPE values")

## Combine the plots 
fig_c1 = ggarrange(fig_c1_a, fig_c1_b, fig_c1_c, 
                  common.legend = TRUE,
                  nrow = 1, ncol = 3,
                  align = "h")
ggsave("figures/FigureC1_SimulationStudy_Varyingq.pdf", fig_c1, 
       device = "pdf", width = 13.04, height = 4.68)

################################################################################
## Load all the case study data 
################################################################################
load("data/Case Study/CaseStudy_EconomicLoss.RData") # results_economic_loss
load("data/Case Study/CaseStudy_SystemBlackout.RData") # results_system_blackout
load("data/Case Study/CaseStudy_RampLambda.RData") # results_ramp_lambda
load("data/Case Study/CaseStudy_RampGamma.RData") # results_ramp_gamma
load_dat = read.csv("data/GEFCom2012/Load_history.csv")
temp_dat = read.csv("data/GEFCom2012/Temperature_history.csv")
# The data below is properly commented in Paper_Code_Data.R.
select = dplyr::select
load_dat = load_dat %>% mutate(zone_id = zone_id %>% as.factor,
                               year = year,
                               month = month %>% as.factor,
                               day = (mod(day - 1, 7) + 1) %>% as.factor)
temp_dat %<>% mutate(station_id = station_id %>% as.factor,
                     year = year,
                     month = month %>% as.factor)
load_dat[, 5:28] = Vectorize(eeptools::decomma)(load_dat[, 5:28])
temp_dat[, 5:28] = Vectorize(eeptools::decomma)(temp_dat[, 5:28])
load_dat = load_dat %>% as_tibble %>%  filter(year >= 2005 & year <= 2007)
temp_dat = temp_dat %>% as_tibble %>% filter(year >= 2005 & year <= 2007)
load_dat %<>% pivot_longer(h1:h24, values_to = "Load", names_to = "hour", names_prefix = "h")
yt = load_dat %>%
  split(load_dat$zone_id) %>%
  flatten %>%
  as.data.frame %>%
  select_if(stringr::str_detect(names(.), "Load")) %>%
  rowSums
load_dat = load_dat %>%
  split(load_dat$zone_id) %>%
  flatten %>%
  as.data.frame %>%
  select(year:Load)
load_dat$Load = yt
temp_dat %<>% pivot_longer(h1:h24, values_to = "Temperature", names_to = "hour", names_prefix = "h")
T = temp_dat %>%
  split(temp_dat$station_id) %>%
  flatten %>%
  as.data.frame %>%
  select_if(stringr::str_detect(names(.), "Temperature")) %>%
  rowMeans
load_dat$Temperature = T
load_dat$Trend = 1:nrow(load_dat)
names(load_dat) = c("Year", "Month", "Weekday", "Hour", "Load", "Temperature", "Trend")
load_dat$Hour %<>% as.factor()
load_dat$Load = load_dat$Load/1e6
load_dat %<>% na.omit()
load_dat_training = load_dat %>% as_tibble %>% filter(Year <= 2006)
load_dat_test = load_dat %>% as_tibble %>% filter(Year == 2007)

################################################################################
## Process the above results into a format appropriate for plotting
################################################################################

## results_economic_loss
all_results_economic_loss = process_all_results(results_economic_loss$results, results_economic_loss$grid_vals) %>%
  tibble %>%
  group_by(Method, mu, cv, q) 
all_results_economic_loss$p = 1 - all_results_economic_loss$p 
all_results_economic_loss %<>% mutate(mu = as.factor(mu), cv = as.factor(cv), Method = as.factor(Method))
levels(all_results_economic_loss$Method)[1] = "Jiao"
all_results_economic_loss %<>% filter(Method %in% c("Jiao", "Bisquare", "Bisquare DD", "Huber", "Huber DD", "New", "LS"))
levels(all_results_economic_loss$cv) = c("1/6", "1/5", "1/4")

## results_system_blackout
all_results_system_blackout = process_all_results(results_system_blackout$results, results_system_blackout$grid_vals) %>%
  tibble %>%
  group_by(Method, mu, cv, q) 
all_results_system_blackout$p = 1 - all_results_system_blackout$p 
all_results_system_blackout %<>% mutate(mu = as.factor(mu), cv = as.factor(cv), Method = as.factor(Method))
levels(all_results_system_blackout$Method)[1] = "Jiao"
all_results_system_blackout %<>% filter(Method %in% c("Jiao", "Bisquare", "Bisquare DD", "Huber", "Huber DD", "New", "LS"))
levels(all_results_system_blackout$cv) = c("-1/6", "-1/5", "-1/4")

## results_ramp_lambda 
all_results_ramp_lambda  = process_all_results(results_ramp_lambda$results, results_ramp_lambda$grid_vals) %>%
  tibble %>%
  group_by(Method, L, lambda, q)
all_results_ramp_lambda$p = 1 - all_results_ramp_lambda$p
all_results_ramp_lambda %<>% mutate(L = as.factor(L), lambda = as.factor(lambda), Method = as.factor(Method))
levels(all_results_ramp_lambda$Method)[1] = "Jiao"
all_results_ramp_lambda %<>% filter(Method %in% c("Jiao", "Bisquare", "Bisquare DD", "Huber", "Huber DD", "New", "LS"))

## results_ramp_gamma 
all_results_ramp_gamma  = process_all_results(results_ramp_gamma$results, results_ramp_gamma$grid_vals) %>%
  tibble %>%
  group_by(Method, L, gamma, q)
all_results_ramp_gamma$p = 1 - all_results_ramp_gamma$p
all_results_ramp_gamma %<>% mutate(L = as.factor(L), gamma = as.factor(gamma), Method = as.factor(Method))
levels(all_results_ramp_gamma$Method)[1] = "Jiao"
all_results_ramp_gamma %<>% filter(Method %in% c("Jiao", "Bisquare", "Bisquare DD", "Huber", "Huber DD", "New", "LS"))
levels(all_results_ramp_gamma$L) = c("L: 100", "L: 200", "L: 300")
levels(all_results_ramp_gamma$gamma) = c(expression(paste(gamma, ": 2")), expression(paste(gamma, ": 3")), expression(paste(gamma, ": 4")))

################################################################################
## Figure 2: Types of attacks.
################################################################################
mu = 40
sigma = 7
L = 50
lambda = 0.03
p = 0.3
set.seed(25611)
clean_data = load_dat 
random_attack_data = random_attack(mu, sigma, clean_data$Load, p)
random_attack_data = clean_data %>% mutate(Load = random_attack_data)
ramp_attack_data = ramp_attack(clean_data$Trend, clean_data$Load, p, lambda, L)
ramp_attack_data = clean_data %>% mutate(Load = ramp_attack_data)

fig_2_a = ggplot((clean_data %>% filter(Year == 2005 & Month == 10))[1:168,], aes(Trend, Load)) %>% plot_aes() + 
  labs(y = "Load (MW)") + 
  geom_line() + 
  ylim(c(0.5, 3.0)) + 
  ggtitle("(a): Clean data")
fig_2_b = ggplot((random_attack_data %>% filter(Year == 2005 & Month == 10))[1:168,], aes(Trend, Load)) %>% plot_aes() + 
  labs(y = "Load (MW)") + 
  geom_line() + 
  ylim(c(0.5, 3.0)) + 
  ggtitle("(b): Random attack")
fig_2_c = ggplot((ramp_attack_data %>% filter(Year == 2005 & Month == 10))[1:168,], aes(Trend, Load)) %>% plot_aes() + 
  labs(y = "Load (MW)") + 
  geom_line() + 
  ylim(c(0.5, 3.0)) + 
  ggtitle("(c): Ramp attack")
fig_2 = ggarrange(fig_2_a, fig_2_b, fig_2_c, 
                   common.legend = TRUE,
                   nrow = 1, ncol = 3,
                   align = "h")
ggsave("figures/Figure2_CaseStudy_Attacks.pdf", fig_2, 
       device = "pdf", width = 13.04, height = 4.68)

################################################################################
## Figure 3: Economic loss results.
################################################################################
fig_3_plots = list()
plot_idx = 1
titles = c(expression("(a): "*mu*" = 50, cv = 1/6"),
           expression("(b): "*mu*" = 50, cv = 1/5"),
           expression("(c): "*mu*" = 50, cv = 1/4"),
           expression("(d): "*mu*" = 100, cv = 1/6"),
           expression("(e): "*mu*" = 100, cv = 1/5"),
           expression("(f): "*mu*" = 100, cv = 1/4"),
           expression("(g): "*mu*" = 150, cv = 1/6"),
           expression("(h): "*mu*" = 150, cv = 1/5"),
           expression("(i): "*mu*" = 150, cv = 1/4"))
for (i in 1:length(levels(all_results_economic_loss$mu))) {
  for (j in 1:length(levels(all_results_economic_loss$cv))) {
    plot_data = all_results_economic_loss %>% 
      filter(mu == levels(all_results_economic_loss$mu)[i],
             cv == levels(all_results_economic_loss$cv)[j])
    fig_3_plots[[plot_idx]] = ggplot(plot_data, aes(p, mape, color = Method, linetype = Method)) %>% plot_aes() +
      geom_line(size = 2) +
      scale_y_continuous(labels = scales::percent, limits = c(min(plot_data$mape), max(plot_data$mape))) +
      labs(y = "MAPE") +
      ggtitle(titles[plot_idx]) +
      scale_linetype_manual(values = c(Jiao = 1, Bisquare = 1, `Bisquare DD` = 4, 
                                       Huber = 1, `Huber DD` = 4, LS = 1, New = 1)) 
    plot_idx = plot_idx + 1
  }
}
fig_3 = ggarrange(plotlist = fig_3_plots, 
                  common.legend = TRUE,
                  nrow = 3, ncol = 3,
                  align = "hv")
ggsave("figures/Figure3_CaseStudy_EconomicLoss.pdf", fig_3, 
       device = "pdf", width = 16.04, height = 12.68)

################################################################################
## Figure 4: System blackout results.
################################################################################
fig_4_plots = list()
plot_idx = 1
titles = c(expression("(a): "*mu*" = -20, cv = -1/6"),
           expression("(b): "*mu*" = -20, cv = -1/5"),
           expression("(c): "*mu*" = -20, cv = -1/4"),
           expression("(d): "*mu*" = -40, cv = -1/6"),
           expression("(e): "*mu*" = -40, cv = -1/5"),
           expression("(f): "*mu*" = -40, cv = -1/4"),
           expression("(g): "*mu*" = -60, cv = -1/6"),
           expression("(h): "*mu*" = -60, cv = -1/5"),
           expression("(i): "*mu*" = -60, cv = -1/4"))
for (i in 1:length(levels(all_results_system_blackout$mu))) {
  for (j in 1:length(levels(all_results_system_blackout$cv))) {
    plot_data = all_results_system_blackout %>% 
      filter(mu == levels(all_results_system_blackout$mu)[i],
             cv == levels(all_results_system_blackout$cv)[j])
    fig_4_plots[[plot_idx]] = ggplot(plot_data, aes(p, mape, color = Method, linetype = Method)) %>% plot_aes() +
      geom_line(size = 2) +
      scale_y_continuous(labels = scales::percent, limits = c(min(plot_data$mape), max(plot_data$mape))) +
      labs(y = "MAPE") +
      ggtitle(titles[plot_idx]) +
      scale_linetype_manual(values = c(Jiao = 1, Bisquare = 1, `Bisquare DD` = 4, 
                                       Huber = 1, `Huber DD` = 4, LS = 1, New = 1))  
    plot_idx = plot_idx + 1
  }
}
fig_4 = ggarrange(plotlist = fig_4_plots, 
                  common.legend = TRUE,
                  nrow = 3, ncol = 3,
                  align = "hv")
ggsave("figures/Figure4_CaseStudy_SystemBlackout.pdf", fig_4, 
       device = "pdf", width = 16.04, height = 12.68)

################################################################################
## Figure 5: Ramp attack results; lambda formulation.
################################################################################
fig_5_plots = list()
plot_idx = 1
titles = c(expression("(a): "*L*" = 40, "*lambda*" = 0.05"),
           expression("(b): "*L*" = 40, "*lambda*" = 0.10"),
           expression("(c): "*L*" = 40, "*lambda*" = 0.15"),
           expression("(d): "*L*" = 60, "*lambda*" = 0.05"),
           expression("(e): "*L*" = 60, "*lambda*" = 0.10"),
           expression("(f): "*L*" = 60, "*lambda*" = 0.15"),
           expression("(g): "*L*" = 80, "*lambda*" = 0.05"),
           expression("(h): "*L*" = 80, "*lambda*" = 0.10"),
           expression("(i): "*L*" = 80, "*lambda*" = 0.15"))
for (i in 1:length(levels(all_results_ramp_lambda$L))) {
  for (j in 1:length(levels(all_results_ramp_lambda$lambda))) {
    plot_data = all_results_ramp_lambda %>% 
      filter(L == levels(all_results_ramp_lambda$L)[i],
             lambda == levels(all_results_ramp_lambda$lambda)[j])
    fig_5_plots[[plot_idx]] = ggplot(plot_data, aes(p, mape, color = Method, linetype = Method)) %>% plot_aes() +
      geom_line(size = 2) +
      scale_y_continuous(labels = scales::percent, limits = c(min(plot_data$mape), max(plot_data$mape))) +
      labs(y = "MAPE") +
      ggtitle(titles[plot_idx]) +
      scale_linetype_manual(values = c(Jiao = 1, Bisquare = 1, `Bisquare DD` = 4, 
                                       Huber = 1, `Huber DD` = 4, LS = 1, New = 1))  
    plot_idx = plot_idx + 1
  }
}
fig_5 = ggarrange(plotlist = fig_5_plots, 
                  common.legend = TRUE,
                  nrow = 3, ncol = 3,
                  align = "hv")
ggsave("figures/Figure5_CaseStudy_RampLambda.pdf", fig_5, 
       device = "pdf", width = 16.04, height = 12.68)

################################################################################
## Figure 6: Ramp attack results; gamma formulation.
################################################################################
fig_6_plots = list()
plot_idx = 1
titles = c(expression("(a): "*L*" = 100, "*gamma*" = 2"),
           expression("(b): "*L*" = 100, "*gamma*" = 3"),
           expression("(c): "*L*" = 100, "*gamma*" = 4"),
           expression("(d): "*L*" = 200, "*gamma*" = 2"),
           expression("(e): "*L*" = 200, "*gamma*" = 3"),
           expression("(f): "*L*" = 200, "*gamma*" = 4"),
           expression("(g): "*L*" = 300, "*gamma*" = 2"),
           expression("(h): "*L*" = 300, "*gamma*" = 3"),
           expression("(i): "*L*" = 300, "*gamma*" = 4"))
for (i in 1:length(levels(all_results_ramp_gamma$L))) {
  for (j in 1:length(levels(all_results_ramp_gamma$gamma))) {
    plot_data = all_results_ramp_gamma %>% 
      filter(L == levels(all_results_ramp_gamma$L)[i],
             gamma == levels(all_results_ramp_gamma$gamma)[j])
    fig_6_plots[[plot_idx]] = ggplot(plot_data, aes(p, mape, color = Method, linetype = Method)) %>% plot_aes() +
      geom_line(size = 2) +
      scale_y_continuous(labels = scales::percent, limits = c(min(plot_data$mape), max(plot_data$mape))) +
      labs(y = "MAPE") +
      ggtitle(titles[plot_idx]) +
      scale_linetype_manual(values = c(Jiao = 1, Bisquare = 1, `Bisquare DD` = 4, 
                                       Huber = 1, `Huber DD` = 4, LS = 1, New = 1))  
    plot_idx = plot_idx + 1
  }
}
fig_6 = ggarrange(plotlist = fig_6_plots, 
                  common.legend = TRUE,
                  nrow = 3, ncol = 3,
                  align = "hv")
ggsave("figures/Figure6_CaseStudy_RampGamma.pdf", fig_6, 
       device = "pdf", width = 16.04, height = 12.68)