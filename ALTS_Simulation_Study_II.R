################################################################################
## Simulation Study II: What is the best value of q? Assessing cross-validation results.
source("ALTS.R")
source("Plotting.R")
library(MLmetrics)
library(pracma)

oldw <- getOption("warn")
options(warn = -1)

## Define the simulation function
sim_fnc_2 = function(beta, X, p, sigma_1, sigma_2) {
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
  ## Fit the results 
  fitControl = trainControl(method = "cv", number = 5, search = "random")
  altsModel = train(model_formula, data = data, metric = "MAE", method = lpALTS, tuneLength = 10, trControl = fitControl)
  ## Return results 
  return(altsModel)
}

## Define function for processing simulation results 
`%notin%` <- Negate(`%in%`)
process_results = function(results, beta, X) {
  p_sim_new = c()
  sigma_sim_new = c()
  mape_sim_new = c()
  q_sim_new = c()
  mu = beta[1] + X %*% beta[2:length(beta)]
  for (i in 1:length(results)) {
    p_sim_new = c(p_sim_new, results[[i]]$finalModel$p[length(results[[i]]$finalModel$p)])
    sigma_sim_new = c(sigma_sim_new, results[[i]]$finalModel$sigma[length(results[[i]]$finalModel$sigma)])
    beta_new = results[[i]]$finalModel$beta[names(results[[i]]$finalModel$beta) %notin% "weights"]
    mu_new = beta_new[1] + X %*% beta_new[2:length(beta_new)]
    mape_sim_new = c(mape_sim_new, MAPE(mu_new, mu))
    q_sim_new = c(q_sim_new, results[[i]]$finalModel$tuneValue$q)
  }
  p_mean_new = mean(p_sim_new, na.rm = TRUE)
  sigma_mean_new = mean(sigma_sim_new, na.rm = TRUE)
  mape_mean_new = mean(mape_sim_new, na.rm = TRUE)
  q_sim_new = mean(q_sim_new, na.rm = TRUE)
  new_results = list(p = p_mean_new, sigma = sigma_mean_new, rmape = mape_mean_new, q = q_sim_new)
  return(list(New = new_results))
}

## Define the simulation function for a single p
main_sim_fnc = function(num_sims, beta, X, p, sigma_1, sigma_2) {
  results = list()
  i = 1
  while (i <= num_sims) {
    tryit = try(sim_fnc_2(beta, X, p, sigma_1, sigma_2), silent = FALSE)
    if (inherits(tryit, "try-error")) {
      print("Error occurred. Restarting loop.")
    } else {
      results[[i]] = tryit
      i = i + 1
    }
  }
  sim_results = process_results(results, beta, X)
  return(sim_results)
}

## Define the complete simulation function for a range of p
complete_sim_fnc = function(num_sims, beta, X, p, sigma_1, sigma_2) {
  new_p = matrix(data = NA, nrow = 1, ncol = length(p))
  new_sigma = matrix(data = NA, nrow = 1, ncol = length(p))
  new_rmape = matrix(data = NA, nrow = 1, ncol = length(p))
  new_q = matrix(data = NA, nrow = 1, ncol = length(p))
  for (j in 1:length(p)) {
    results = main_sim_fnc(num_sims, beta, X, p[j], sigma_1, sigma_2)
    new_p[1, j] = results$New$p
    new_sigma[1, j] = results$New$sigma
    new_rmape[1, j] = results$New$rmape
    new_q[1, j] = results$New$q
    print(j)
  }
  new_results = list(p = new_p, sigma = new_sigma, rmape = new_rmape, q = new_q)
  return(new_results)
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
sigma_1 = 0.1
sigma_2 = 10 

## Results 
set.seed(1003991)
results = complete_sim_fnc(num_sims, beta, X, p, sigma_1, sigma_2)

options(warn = oldw)

## Collect results 
new_p = results$p %>% t
new_sigma = results$sigma %>% t
new_rmape = results$rmape %>% t
new_q = results$q %>% t

library(tidyverse)

new_p %<>% as.data.frame %>% mutate(p = p) 
names(new_p)[1] = "phat"
new_sigma %<>% as.data.frame %>% mutate(p = p)
names(new_sigma)[1] = "sigmahat"
new_rmape  %<>% as.data.frame %>% mutate(p = p) 
names(new_rmape)[1] = "mape"
new_q  %<>% as.data.frame %>% mutate(p = p) 
names(new_q)[1] = "qhat"

## Plot 
library(ggplot2)
library(lemon)
library(scales)
plt1 = ggplot(new_p, aes(p, phat)) %>% plot_aes() + geom_line(size = 2) +
  labs(y = expression(hat(p))) + 
  geom_abline(slope = 1, size = 1, color = "red", linetype = "dashed") +
  xlim(c(0.5, 1)) + ylim(c(0.5, 1))

plt2 = ggplot(new_sigma, aes(p, sigmahat)) %>% plot_aes() + geom_line(size = 2) +
  labs(y = expression(hat(sigma))) + 
  geom_hline(yintercept = sigma_1, color = "red", linetype = "dashed") +
  ylim(c(0.09, 0.2))

plt3 = ggplot(new_rmape, aes(p, mape)) %>% plot_aes() + geom_line(size = 2) +
  labs(y = "MAPE") 
  ylim(c(0, 0.025))
  
plt4 = ggplot(new_q, aes(p, qhat)) %>% plot_aes() + geom_line(size = 2) + 
  labs(y = expression(hat(q)))

plts = ggarrange(plt1, plt2, plt3, plt4, nrow = 2, ncol = 2)
plts

## Comparing to Study I
load("crossvalidation_results.RData")
results_II = results 
load("varyingq_results.RData")
results_I = results

p = c(0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)
q = c(1.0, 1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14)
new_p = results_I$New$p %>% t
new_sigma = results_I$New$sigma %>% t
new_rmape = results_I$New$rmape %>% t
bacher_p = results_I$Bacher$p %>% t
bacher_sigma = results_I$Bacher$sigma %>% t
bacher_rmape = results_I$Bacher$rmape %>% t
new_p %<>% as.data.frame %>% mutate(p = p) %>% pivot_longer(cols = 1:length(q), names_to = "q", values_to = "phat")
new_sigma %<>% as.data.frame %>% mutate(p = p) %>% pivot_longer(cols = 1:length(q), names_to = "q", values_to = "sigmahat")
new_rmape  %<>% as.data.frame %>% mutate(p = p) %>% pivot_longer(cols = 1:length(q), names_to = "q", values_to = "mape")
bacher_p %<>% as.data.frame %>% mutate(p = p)
names(bacher_p)[1] = "phat"
bacher_sigma %<>% as.data.frame %>% mutate(p = p)
names(bacher_sigma)[1] = "sigmahat"
bacher_rmape %<>% as.data.frame %>% mutate(p = p)
names(bacher_rmape)[1] = "mape"


plt1 = ggplot(new_p, aes(p, phat, color = q)) %>% plot_aes() + geom_line(size = 2, alpha = 0.3) +
  scale_color_discrete(labels = q) + labs(y = expression(hat(p))) + 
  geom_abline(slope = 1, size = 2) +
  geom_line(data = bacher_p, aes(p, phat), color = "red", linetype = "dashed", size = 1, alpha = 0.3) +
  xlim(c(0.5, 1)) + ylim(c(0.5, 1))

plt2 = ggplot(new_sigma, aes(p, sigmahat, color = q)) %>% plot_aes() + geom_line(size = 2, alpha = 0.3) +
  scale_color_discrete(labels = q) + labs(y = expression(hat(sigma))) + 
  geom_hline(yintercept = sigma_1) +
  geom_line(data = bacher_sigma, aes(p, sigmahat), color = "red", linetype = "dashed", size = 1, alpha = 0.3) + ylim(c(0.09, 0.2))

plt3 = ggplot(new_rmape, aes(p, mape, color = q)) %>% plot_aes() + geom_line(size = 2, alpha = 0.3) +
  scale_color_discrete(labels = q) + labs(y = "MAPE") + scale_y_continuous(labels=scales::percent) +
  geom_line(data = bacher_rmape, aes(p, mape), color = "red", linetype = "dashed", size = 2, alpha = 0.3) + ylim(c(0, 0.025))

new_p = results_II$p %>% t
new_sigma = results_II$sigma %>% t
new_rmape = results_II$rmape %>% t
new_q = results_II$q %>% t
new_p %<>% as.data.frame %>% mutate(p = p) 
names(new_p)[1] = "phat"
new_sigma %<>% as.data.frame %>% mutate(p = p)
names(new_sigma)[1] = "sigmahat"
new_rmape  %<>% as.data.frame %>% mutate(p = p) 
names(new_rmape)[1] = "mape"
new_q  %<>% as.data.frame %>% mutate(p = p) 
names(new_q)[1] = "qhat"

plt1 = plt1 + geom_line(data = new_p, aes(p, phat), color = "blue", size = 2)
plt2 = plt2 + geom_line(data = new_sigma, aes(p, sigmahat), color = "blue", size = 2)
plt3 = plt3 + geom_line(data = new_rmape, aes(p, mape), color = "blue", size = 2)

plts = ggarrange(plt1, plt2, plt3, common.legend = TRUE, nrow = 1, ncol = 3, labels = c("(a)", "(b)", "(c)"))
ggsave("Figure2_CrossValidation.pdf", plts, device = "pdf", width = 13.04, height = 4.68)