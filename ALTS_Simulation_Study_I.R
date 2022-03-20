################################################################################
## Simulation Study I: Adding some noise.
source("ALTS.R")
source("Plotting.R")
library(MLmetrics)
library(pracma)

oldw <- getOption("warn")
options(warn = -1)

## Define the simulation function
sim_fnc_1 = function(beta, X, p, sigma_1, sigma_2, q, include_bacher = TRUE) {
  ## Add noise to the data 
  nobs = nrow(X)
  mu = beta[1] + X %*% beta[2:length(beta)]
  h = floor(p * nobs)
  noise_good = rnorm(h, 0, sigma_1)
  noise_outlier = rnorm(nobs - h, 0, sigma_2)
  # noise_outlier = rep(60, nobs - h)
  noise = c(noise_good, noise_outlier)
  y = mu + noise
  data = as.data.frame(X)
  data$y = y 
  model_formula = y ~ .
  ## Bacher results
  if (include_bacher) {
    bacher_results = ALTS_Bacher(model_formula, data)
  } else {
    bacher_results = NA
  }
  ## New results
  new_results = ALTS(model_formula, data, q)
  ## Return results 
  return(list(Bacher = bacher_results, New = new_results))
}

## Define function for processing simulation results
`%notin%` <- Negate(`%in%`)
process_results = function(results, beta, X, include_bacher = TRUE) {
  if (include_bacher) {
    p_sim_bacher = c()
    sigma_sim_bacher = c()
    mape_sim_bacher = c()
  }
  p_sim_new = c()
  sigma_sim_new = c()
  mape_sim_new = c()
  mu = beta[1] + X %*% beta[2:length(beta)]
  for (i in 1:length(results)) {
    if (include_bacher) {
      p_sim_bacher = c(p_sim_bacher, results[[i]]$Bacher$p[length(results[[i]]$Bacher$p)])
      sigma_sim_bacher = c(sigma_sim_bacher, results[[i]]$Bacher$sigma[length(results[[i]]$Bacher$sigma)])
      beta_bacher = results[[i]]$Bacher$beta[names(results[[i]]$Bacher$beta) %notin% "weights"]
      mu_bacher = beta_bacher[1] + X %*% beta_bacher[2:length(beta_bacher)] 
      mape_sim_bacher = c(mape_sim_bacher, MAPE(mu_bacher, mu))
    }
    p_sim_new = c(p_sim_new, results[[i]]$New$p[length(results[[i]]$New$p)])
    sigma_sim_new = c(sigma_sim_new, results[[i]]$New$sigma[length(results[[i]]$New$sigma)])
    beta_new = results[[i]]$New$beta[names(results[[i]]$New$beta) %notin% "weights"]
    mu_new = beta_new[1] + X %*% beta_new[2:length(beta_new)]
    mape_sim_new = c(mape_sim_new, MAPE(mu_new, mu))
  }
  if (include_bacher) {
    p_mean_bacher = mean(p_sim_bacher, na.rm = TRUE)
    sigma_mean_bacher = mean(sigma_sim_bacher, na.rm = TRUE)
    mape_mean_bacher = mean(mape_sim_bacher, na.rm = TRUE)
    bacher_results = list(p = p_mean_bacher, sigma = sigma_mean_bacher, rmape = mape_mean_bacher)
  } else {
    bacher_results = NA
  }
  p_mean_new = mean(p_sim_new, na.rm = TRUE)
  sigma_mean_new = mean(sigma_sim_new, na.rm = TRUE)
  mape_mean_new = mean(mape_sim_new, na.rm = TRUE)
  new_results = list(p = p_mean_new, sigma = sigma_mean_new, rmape = mape_mean_new)
  return(list(Bacher = bacher_results, New = new_results))
}

## Define the total simulation function for a single (p, q)
main_sim_fnc = function(num_sims, beta, X, p, sigma_1, sigma_2, q, include_bacher = TRUE) {
  results = list()
  i = 1
  while (i <= num_sims) {
    tryit = try(sim_fnc_1(beta, X, p, sigma_1, sigma_2, q, include_bacher = include_bacher), silent = FALSE)
    if (inherits(tryit, "try-error")) {
      print("Error occurred. Restarting loop.")
    } else {
      results[[i]] = tryit
      i = i + 1
    }
  }
  sim_results = process_results(results, beta, X, include_bacher = include_bacher)
  return(sim_results)
}

## Define the complete simulation function for a range of (p, q)
complete_sim_fnc = function(num_sims, beta, X, p, sigma_1, sigma_2, q) {
  bacher_p = matrix(data = NA, nrow = 1, ncol = length(p))
  bacher_sigma = matrix(data = NA, nrow = 1, ncol = length(p))
  bacher_rmape = matrix(data = NA, nrow = 1, ncol = length(p))
  new_p = matrix(data = NA, nrow = length(q), ncol = length(p))
  new_sigma = matrix(data = NA, nrow = length(q), ncol = length(p))
  new_rmape = matrix(data = NA, nrow = length(q), ncol = length(p))
  for (j in 1:length(p)) {
    results = main_sim_fnc(num_sims, beta, X, p[j], sigma_1, sigma_2, q[1], include_bacher = TRUE)
    bacher_p[1, j] = results$Bacher$p
    bacher_sigma[1, j] = results$Bacher$sigma
    bacher_rmape[1, j] = results$Bacher$rmape
    new_p[1, j] = results$New$p
    new_sigma[1, j] = results$New$sigma
    new_rmape[1, j] = results$New$rmape
    print(c(1, j))
  }
  for (i in 2:length(q)) {
    for (j in 1:length(p)) {
      results = main_sim_fnc(num_sims, beta, X, p[j], sigma_1, sigma_2, q[i], include_bacher = FALSE)
      new_p[i, j] = results$New$p
      new_sigma[i, j] = results$New$sigma
      new_rmape[i, j] = results$New$rmape
      print(c(i, j))
    }
  }
  bacher_results = list(p = bacher_p, sigma = bacher_sigma, rmape = bacher_rmape)
  new_results = list(p = new_p, sigma = new_sigma, rmape = new_rmape)
  return(list(Bacher = bacher_results, New = new_results))
}

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
num_sims = 100
p = c(0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
q = c(1.0, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14)
pq_grid_p = meshgrid(p, q)$X 
pq_grid_q = meshgrid(p, q)$Y
sigma_1 = 0.1
sigma_2 = 10 

## Results 
set.seed(1003991)
results = complete_sim_fnc(num_sims, beta, X, p, sigma_1, sigma_2, q)

options(warn = oldw)

## Collect results 
new_p = results$New$p %>% t
new_sigma = results$New$sigma %>% t
new_rmape = results$New$rmape %>% t
bacher_p = results$Bacher$p %>% t
bacher_sigma = results$Bacher$sigma %>% t
bacher_rmape = results$Bacher$rmape %>% t

library(tidyverse)

new_p %<>% as.data.frame %>% mutate(p = p) %>% pivot_longer(cols = 1:length(q), names_to = "q", values_to = "phat")
new_sigma %<>% as.data.frame %>% mutate(p = p) %>% pivot_longer(cols = 1:length(q), names_to = "q", values_to = "sigmahat")
new_rmape  %<>% as.data.frame %>% mutate(p = p) %>% pivot_longer(cols = 1:length(q), names_to = "q", values_to = "mape")
bacher_p %<>% as.data.frame %>% mutate(p = p)
names(bacher_p)[1] = "phat"
bacher_sigma %<>% as.data.frame %>% mutate(p = p)
names(bacher_sigma)[1] = "sigmahat"
bacher_rmape %<>% as.data.frame %>% mutate(p = p)
names(bacher_rmape)[1] = "mape"

## Plot 
library(ggplot2)
library(lemon)
library(scales)
plt1 = ggplot(new_p, aes(p, phat, color = q)) %>% plot_aes() + geom_line(size = 2) +
  scale_color_discrete(labels = q) + labs(y = expression(hat(p))) + 
  geom_abline(slope = 1, size = 2) +
  geom_line(data = bacher_p, aes(p, phat), color = "red", linetype = "dashed", size = 1) +
  xlim(c(0.49, 1)) + ylim(c(0.49, 1))

plt2 = ggplot(new_sigma, aes(p, sigmahat, color = q)) %>% plot_aes() + geom_line(size = 2) +
  scale_color_discrete(labels = q) + labs(y = expression(hat(sigma))) + 
  geom_hline(yintercept = sigma_1, size = 2) +
  geom_line(data = bacher_sigma, aes(p, sigmahat), color = "red", linetype = "dashed", size = 1) + ylim(c(0.09, 0.25))

plt3 = ggplot(new_rmape, aes(p, mape, color = q)) %>% plot_aes() + geom_line(size = 2) +
  scale_color_discrete(labels = q) + labs(y = "MAPE") + 
  geom_line(data = bacher_rmape, aes(p, mape), color = "red", linetype = "dashed", size = 2) + ylim(c(0, 0.025)) +
  scale_y_continuous(labels=scales::percent)

plts = ggarrange(plt1, plt2, plt3, common.legend = TRUE, nrow = 1, ncol = 3, labels = c("(a)", "(b)", "(c)"))
ggsave("Figure1_VaryingQ.pdf", plts, device = "pdf", width  = 13.04, height = 4.68)


