################################################################################
###### Load the data
################################################################################
library(tidyverse)
library(magrittr)
library(eeptools)
library(dplyr)
library(ggpubr)
library(MLmetrics)
library(ggplot2)
library(quantreg)
source("Cyberattack_Functions.R")
source("Plotting.R")
source("ALTS.R")
load_dat = read.csv("data/Load_history.csv")
temp_dat = read.csv("data/Temperature_history.csv")
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

# Fix hour 
load_dat$Hour %<>% as.factor()

# Change names 
names(load_dat) = c("Year", "Month", "Weekday", "Hour", "Load", "Temperature", "Trend")

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
###### Looking at the data 
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

clean_plot = ggplot((clean_data %>% filter(Year == 2005 & Month == 10))[1:168,], aes(Trend, Load)) %>% plot_aes() + 
  labs(y = "Load (MW)") + geom_line() + ylim(c(0.5, 3.0)) + ggtitle("(a): Clean data")
random_plot = ggplot((random_attack_data %>% filter(Year == 2005 & Month == 10))[1:168,], aes(Trend, Load)) %>% plot_aes() + 
  labs(y = "Load (MW)") + geom_line() + ylim(c(0.5, 3.0)) + ggtitle("(b): Random attack")
ramp_plot = ggplot((ramp_attack_data %>% filter(Year == 2005 & Month == 10))[1:168,], aes(Trend, Load)) %>% plot_aes() + 
  labs(y = "Load (MW)") + geom_line() + ylim(c(0.5, 3.0)) + ggtitle("(c): Ramp attack")
all_plots = ggarrange(clean_plot, random_plot, ramp_plot, nrow = 1, ncol = 3)
ggsave("Figure3_TypesOfAttacks.pdf", all_plots, device = "pdf", width = 13.04, height = 4.68)

################################################################################
###### Economic loss studies 
################################################################################
random_attack_ALTS = function(load_dat_training, load_dat_test, mu, sigma, p, q, model_formula) {
  new_data = load_dat_training %>% mutate(Load = random_attack(mu, sigma, load_dat_training$Load, p))
  results = ALTS(model_formula, new_data, q, use_lad = FALSE)
  p_hat = results$p[length(results$p)]
  print(results$p)
  sigma_hat = results$sigma[length(results$sigma)]
  mape_val = MLmetrics::MAPE(predict(results$model, load_dat_test), load_dat_test$Load)
  new_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  results = ALTS_Bacher(model_formula, new_data, use_lad = FALSE)
  p_hat = results$p[length(results$p)]
  sigma_hat = results$sigma[length(results$sigma)]
  mape_val = MLmetrics::MAPE(predict(results$model, load_dat_test), load_dat_test$Load)
  bacher_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  return(list(New = new_results, Bacher = bacher_results))
}

random_attack_all_simulations = function(num_sims, load_dat_training, load_dat_test, mu, sigma, p, q, model_formula) {
  new_p = c()
  new_sigma = c()
  new_mape = c()
  bacher_p = c()
  bacher_sigma = c()
  bacher_mape = c()
  for (i in 1:num_sims) {
    results = random_attack_ALTS(load_dat_training, load_dat_test, mu, sigma, p, q, model_formula)
    new_p = c(new_p, results$New$p)
    new_sigma = c(new_sigma, results$New$sigma)
    new_mape = c(new_mape, results$New$mape)
    bacher_p = c(bacher_p, results$Bacher$p)
    bacher_sigma = c(bacher_sigma, results$Bacher$sigma)
    bacher_mape = c(bacher_mape, results$Bacher$mape)
  }
  new_p_mean = mean(new_p)
  new_sigma_mean = mean(new_sigma)
  new_mape_mean = mean(new_mape)
  bacher_p_mean = mean(bacher_p)
  bacher_sigma_mean = mean(bacher_sigma)
  bacher_mape_mean = mean(bacher_mape)
  new_results = list(p = new_p_mean, sigma = new_sigma_mean, mape = new_mape_mean)
  bacher_results = list(p = bacher_p_mean, sigma = bacher_sigma_mean, mape = bacher_mape_mean)
  return(list(New = new_results, Bacher = bacher_results))
}

mu = c(50, 100, 150, 200)
sigma = c(1/6, 1/5, 1/4, 1/3) * mu 
cv = sigma/mu
p = seq(0, 0.45, by = 0.05)
q = 1.10 
num_sims = 100
grid_vals = expand.grid(mu = mu, sigma = sigma, p = p)
results = random_attack_all_simulations(3, load_dat_training, load_dat_test, mu[3], sigma[3], p[3], 1.08, model_formula)

