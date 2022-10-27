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
load_dat = read.csv("data/GEFCom2012/Load_history.csv")
temp_dat = read.csv("data/GEFCom2012/Temperature_history.csv")
select = dplyr::select
set.seed(1230001)
load_dat = load_dat %>% mutate(zone_id = zone_id %>% as.factor,year = year,month = month %>% as.factor,day = (mod(day - 1, 7) + 1) %>% as.factor)
temp_dat %<>% mutate(station_id = station_id %>% as.factor,year = year,month = month %>% as.factor)
load_dat[, 5:28] = Vectorize(eeptools::decomma)(load_dat[, 5:28])
temp_dat[, 5:28] = Vectorize(eeptools::decomma)(temp_dat[, 5:28])
load_dat = load_dat %>% as_tibble %>%  filter(year >= 2005 & year <= 2007)
temp_dat = temp_dat %>% as_tibble %>% filter(year >= 2005 & year <= 2007)
load_dat %<>% pivot_longer(h1:h24, values_to = "Load", names_to = "hour", names_prefix = "h")
yt = load_dat %>%split(load_dat$zone_id) %>%flatten %>%as.data.frame %>%select_if(stringr::str_detect(names(.), "Load")) %>%rowSums
load_dat = load_dat %>%split(load_dat$zone_id) %>%flatten %>%as.data.frame %>%select(year:Load)
load_dat$Load = yt
temp_dat %<>% pivot_longer(h1:h24, values_to = "Temperature", names_to = "hour", names_prefix = "h")
T = temp_dat %>%split(temp_dat$station_id) %>%flatten %>%as.data.frame %>%select_if(stringr::str_detect(names(.), "Temperature")) %>%rowMeans
load_dat$Temperature = T
load_dat$Trend = 1:nrow(load_dat)
names(load_dat) = c("Year", "Month", "Weekday", "Hour", "Load", "Temperature", "Trend")
load_dat$Hour %<>% as.factor()
load_dat$Load = load_dat$Load/1e6
load_dat %<>% na.omit()
load_dat_training = load_dat %>% as_tibble %>% filter(Year <= 2006)
load_dat_test = load_dat %>% as_tibble %>% filter(Year == 2007)
model_formula = Load ~ Trend + Month + Weekday + Hour + Weekday:Hour + Temperature + I(Temperature^2) + I(Temperature^3) +Temperature:Month + I(Temperature^2):Month + I(Temperature^3):Month +Temperature:Hour + I(Temperature^2):Hour + I(Temperature^3):Hour

################################################################################
## Economic loss study: Smaller mu, larger sigma 
################################################################################
set.seed(222328391)
mu = c(5, 10, 15)
sigma = c(2, 4, 6) * mu 
cv = sigma/mu
p = seq(0.1, 0.5, by = 0.1)
q = 1.2
num_sims = 5
results_economic_loss_small_mu = random_attack_all_simulations_grid_vals(num_sims, load_dat_training, 
                                                                         load_dat_test,
                                                                mu, cv, p, q, model_formula)
