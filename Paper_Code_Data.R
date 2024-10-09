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
library(conquer)

patch_rq_fit_conquer <- function (x, y, tau = 0.5, kernel = c("Gaussian", "uniform", 
                                                              "parabolic", "triangular"), h = 0, tol = 1e-04, iteMax = 5000, 
                                  ci = FALSE, alpha = 0.05, B = 200) 
{
  if (!requireNamespace("conquer", quietly = TRUE)) 
    stop("method conquer requires package conquer")
  fit = conquer::conquer(x[, -1], y, tau = tau, kernel = kernel, 
                         h = h, tol = tol, iteMax = iteMax, ci = ci, alpha = alpha, 
                         B = 1000)  # <-- here, passing ci as an argument (original version used ci = FALSE)
  coefficients = fit$coeff
  names(coefficients) = dimnames(x)[[2]]
  residuals = fit$residual
  list(coefficients = coefficients, tau = tau, residuals = residuals)
}

# Overriding the function in the package's namespace
assignInNamespace("rq.fit.conquer", patch_rq_fit_conquer, ns = "quantreg")

################################################################################
## Simulation study
################################################################################
## Define the model 
set.seed(99291)
nobs = 2000
beta0 = -1.3
beta1 = 2.0
beta2 = 1.7
beta3 = -3.0
x1 = runif(nobs, -1, 1)
x2 = rnorm(nobs, 0, 1)
x3 = runif(nobs, 0, 1)
beta = c(beta0, beta1, beta2, beta3)
X = cbind(x1, x2, x3)

## Parameters 
sigma_1 = 0.1
sigma_2 = 1.3

## Results: fixed q 
set.seed(1003991)
num_sims_1 = 100
p_1 = seq(0.5, 1.0, by = 0.05)
q_1 = 1.35
results_fixed_q = complete_sim_fnc(num_sims_1, beta, X, p_1, sigma_1, sigma_2, q = q_1)
save(results_fixed_q, file = "data/Simulation/Fixedq.RData")

## Results: varying q
set.seed(100399391)
num_sims_2 = 100
p_2 = seq(0.5, 1.0, by = 0.05)
q_2 = seq(1.0, 1.5, by = 0.07)
results_varying_q = complete_sim_fnc_vary_q(num_sims_2, beta, X, p_2, sigma_1, sigma_2, q_2)
save(results_varying_q, file = "data/Simulation/Varyingq.RData")

################################################################################
## Read in the case study data
################################################################################
load_dat = read.csv("data/GEFCom2012/Load_history.csv")
temp_dat = read.csv("data/GEFCom2012/Temperature_history.csv")
select = dplyr::select
set.seed(1230001)

# Fix variable types
load_dat = load_dat %>% mutate(zone_id = zone_id %>% as.factor,
                               year = year,
                               month = month %>% as.factor,
                               day = (mod(day - 1, 7) + 1) %>% as.factor)
temp_dat %<>% mutate(station_id = station_id %>% as.factor,
                     year = year,
                     month = month %>% as.factor)

# Remove commas
load_dat[, 5:28] = Vectorize(eeptools::decomma)(load_dat[, 5:28])
temp_dat[, 5:28] = Vectorize(eeptools::decomma)(temp_dat[, 5:28])

# Filter within a certain time period 
load_dat = load_dat %>% as_tibble %>%  filter(year >= 2005 & year <= 2007)
temp_dat = temp_dat %>% as_tibble %>% filter(year >= 2005 & year <= 2007)

# Melt hours 
load_dat %<>% pivot_longer(h1:h24, values_to = "Load", names_to = "hour", names_prefix = "h")

# Create the response variable: sums of loads across the 24 stations
yt = load_dat %>%
  split(load_dat$zone_id) %>%
  flatten %>%
  as.data.frame %>%
  select_if(stringr::str_detect(names(.), "Load")) %>%
  rowSums

# Fix data 
load_dat = load_dat %>%
  split(load_dat$zone_id) %>%
  flatten %>%
  as.data.frame %>%
  select(year:Load)
load_dat$Load = yt

# Create the temperature variable 
temp_dat %<>% pivot_longer(h1:h24, values_to = "Temperature", names_to = "hour", names_prefix = "h")
T = temp_dat %>%
  split(temp_dat$station_id) %>%
  flatten %>%
  as.data.frame %>%
  select_if(stringr::str_detect(names(.), "Temperature")) %>%
  rowMeans
load_dat$Temperature = T

# Add the linear trend 
load_dat$Trend = 1:nrow(load_dat)

# Change names 
names(load_dat) = c("Year", "Month", "Weekday", "Hour", "Load", "Temperature", "Trend")

# Fix hour 
load_dat$Hour %<>% as.factor()

# Training and testing data
load_dat$Load = load_dat$Load/1e6
load_dat %<>% na.omit()
load_dat_training = load_dat %>% as_tibble %>% filter(Year <= 2006)
load_dat_test = load_dat %>% as_tibble %>% filter(Year == 2007)

# The model 
model_formula = Load ~ Trend + Month + Weekday + Hour + 
  Weekday:Hour + Temperature + I(Temperature^2) + I(Temperature^3) +
  Temperature:Month + I(Temperature^2):Month + I(Temperature^3):Month +
  Temperature:Hour + I(Temperature^2):Hour + I(Temperature^3):Hour

################################################################################
## Economic loss study
################################################################################
set.seed(229591)
mu = c(50, 100, 150)
sigma = c(1/6, 1/5, 1/4) * mu 
cv = sigma/mu
p = seq(0.0, 0.5, by = 0.1)
q = 1.2
num_sims = 5
results_economic_loss = random_attack_all_simulations_grid_vals(num_sims, load_dat_training, load_dat_test,
                                                                mu, cv, p, q, model_formula)
save(results_economic_loss, file = "data/CaseStudy/EconomicLoss.RData")

################################################################################
## Economic loss study: Smaller mu, larger sigma 
################################################################################
set.seed(222328391)
mu = c(5, 10, 15)
sigma = c(2, 4, 6) * mu 
cv = sigma/mu
p = seq(0.0, 0.5, by = 0.1)
q = 1
num_sims = 5
results_economic_loss_small_mu = random_attack_all_simulations_grid_vals(num_sims, load_dat_training, load_dat_test,
                                                                         mu, cv, p, q, model_formula)
save(results_economic_loss_small_mu, file = "data/CaseStudy/EconomicLossSmallMu.RData")

################################################################################
## System blackout study
################################################################################
set.seed(891)
mu = c(-20, -40, -60)
sigma = c(-1/6, -1/5, -1/4) * mu 
cv = sigma/mu
p = seq(0.0, 0.5, by = 0.1)
q = 1.2
num_sims = 5
results_system_blackout = random_attack_all_simulations_grid_vals(num_sims, load_dat_training, load_dat_test,
                                                                mu, cv, p, q, model_formula)
save(results_system_blackout, file = "data/CaseStudy/SystemBlackout.RData")

################################################################################
## Ramp attack: lambda
################################################################################
set.seed(2221)
L = c(40, 60, 80)
lambda = c(0.05, 0.1, 0.15)
p = seq(0.0, 0.5, by = 0.1)
q = 1.2
num_sims = 5
results_ramp_lambda = ramp_attack_all_simulations_grid_vals(num_sims, load_dat_training, load_dat_test,
                                                                  lambda, L, p, q, model_formula)
save(results_ramp_lambda, file = "data/CaseStudy/RampLambda.RData")

################################################################################
## Ramp attack: gamma
################################################################################
set.seed(2666221)
L = c(100, 200, 300)
gamma = c(2, 3, 4)
p = seq(0.0, 0.5, by = 0.1)
q = 1.2
num_sims = 5
results_ramp_gamma = ramp_attack_all_simulations_grid_vals_2(num_sims, load_dat_training, load_dat_test,
                                                            gamma, L, p, q, model_formula)
save(results_ramp_gamma, file = "data/CaseStudy/RampGamma.RData")

