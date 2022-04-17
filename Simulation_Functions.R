################################################################################
## Simulation study functions
################################################################################
## Contents (in order of presentation):
##  1: %notin%
##          An infix operator defining the negation of %in%.
##
##  2. sim_fnc_1(beta, X, p, sigma_1, sigma_2, q = 1.35, include_bacher = TRUE)
##          Perform a single simulation for the simulation study.
##
##  3. process_results(results, beta, X, include_bacher = TRUE)
##          Process the results from a number of simulations from the simulation study.
##
##  4. main_sim_fnc(num_sims, beta, X, p, sigma_1, sigma_2, q, include_bacher = TRUE)
##          For a single `p`, perform the simulation study for `num_sims` simulations.
##
##  5. complete_sim_fnc(num_sims, beta, X, p, sigma_1, sigma_2, q = 1.35)
##          For each `p` in the vector `p`, perform the simulation study for `num_sims` simulations.
##
##  6. complete_sim_fnc_vary_q(num_sims, beta, X, p, sigma_1, sigma_2, q)
##          For each `p` in the vector `p`, and each `q` in the vector `q`, perform the simulation study for `num_sims` simulations.
library(EnvStats)
library(L1pack)
library(dplyr)
library(caret)
library(formula.tools)
library(MLmetrics)
library(pracma)

################################################################################
## %notin%: An infix operator defining the negation of %in%.
################################################################################
`%notin%` <- Negate(`%in%`)

################################################################################
## sim_fnc_1(beta, X, p, sigma_1, sigma_2, q = 1.35, include_bacher = TRUE)
##  Perform a single simulation for the simulation study.
##
## Arguments:
##  beta: Known regression coefficients.
##  X: Design matrix for the regression.
##  p: Proportion of clean data.
##  sigma_1: Standard deviation of the clean data.
##  sigma_2: Standard deviation of the contaminated data.
##  q: The value of q to use in the new ALTS.
##  include_bacher = TRUE: If TRUE, results from Bacher's ALTS should be included.
##
## Outputs:
##  The output is a list containing the following, with each giving the results from `ALTS`:
##    $Bacher: Results from Bacher's ALTS.
##    $New: Results from the new ALTS.
################################################################################
sim_fnc_1 = function(beta, X, p, sigma_1, sigma_2, q = 1.35, include_bacher = TRUE) {
  ## Add noise to the data 
  nobs = nrow(X)
  mu = beta[1] + X %*% beta[2:length(beta)]
  h = floor(p * nobs)
  noise_good = rnorm(h, 0, sigma_1)
  noise_outlier = rnorm(nobs - h, 0, sigma_2)
  noise = c(noise_good, noise_outlier)
  y = mu + noise
  data = as.data.frame(X)
  data$y = as.numeric(y)
  model_formula = y ~ x1 + x2 + x3
  residual_initial = lad(model_formula, data = data)$residuals
  ## Bacher results
  if (include_bacher) {
    bacher_results = ALTS_Bacher(model_formula, data, residual_initial = residual_initial)
  } else {
    bacher_results = NA
  }
  ## New results
  new_results = ALTS(model_formula, data, residual_initial = residual_initial, q = q, maxIters = 20)
  ## Return results 
  return(list(Bacher = bacher_results, New = new_results))
}

################################################################################
## process_results(results, beta, X, include_bacher = TRUE) 
##  Process the results from a number of simulations from the simulation study.
##  The function takes all the individual results and averages over the simulations.
##
## Arguments:
##  results: A list of results of length `num_sims` from `sim_fnc_1`.
##  beta: Known regression coefficients.
##  X: Design matrix for the regression.
##  include_bacher = TRUE: If TRUE, results from Bacher's ALTS should be included.
##
## Outputs:
##  The output is a list containing the following:
##    $Bacher: Results from Bacher's ALTS, further containing:
##      $p: Mean estimates for `p` from Bacher's ALTS.
##      $sigma: Mean estimates for `sigma_1` from Bacher's ALTS.
##      $rmape: Mean MAPE values from Bacher's ALTS.
##    $New: Results from the new ALTS.
##      $q: Value of `q` used from the new ALTS.
##      $p: Mean estimates for `p` from the new ALTS.
##      $sigma: Mean estimates for `sigma_1` from the new ALTS.
##      $rmape: Mean MAPE values from Bacher's ALTS.
################################################################################
process_results = function(results, beta, X, include_bacher = TRUE) {
  if (include_bacher) {
    p_sim_bacher = c()
    sigma_sim_bacher = c()
    mape_sim_bacher = c()
  }
  q_sim_new = c()
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
    q_sim_new = c(q_sim_new, results[[i]]$New$q)
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
  q_mean_new = mean(q_sim_new, na.rm = TRUE)
  p_mean_new = mean(p_sim_new, na.rm = TRUE)
  sigma_mean_new = mean(sigma_sim_new, na.rm = TRUE)
  mape_mean_new = mean(mape_sim_new, na.rm = TRUE)
  new_results = list(p = p_mean_new, sigma = sigma_mean_new, rmape = mape_mean_new, q = q_mean_new)
  return(list(Bacher = bacher_results, New = new_results))
}

################################################################################
## main_sim_fnc(num_sims, beta, X, p, sigma_1, sigma_2, q, include_bacher = TRUE) 
##  For a single `p` and `q`, perform the simulation study for `num_sims` simulations.
##
## Arguments:
##  num_sims: Number of simulations to perform.
##  beta: Known regression coefficients.
##  X: Design matrix for the regression.
##  p: Proportion of clean data.
##  sigma_1: Standard deviation of the clean data.
##  sigma_2: Standard deviation of the contaminated data.
##  q: The value of q to use in the new ALTS.
##  include_bacher = TRUE: If TRUE, results from Bacher's ALTS should be included.
##
## Outputs:
##  The output is a list containing the following, produced from `process_results`:
##    $Bacher: Results from Bacher's ALTS, further containing:
##      $p: Mean estimates for `p` from Bacher's ALTS.
##      $sigma: Mean estimates for `sigma_1` from Bacher's ALTS.
##      $rmape: Mean MAPE values from Bacher's ALTS.
##    $New: Results from the new ALTS.
##      $q: Value of `q` used from the new ALTS.
##      $p: Mean estimates for `p` from the new ALTS.
##      $sigma: Mean estimates for `sigma_1` from the new ALTS.
##      $rmape: Mean MAPE values from Bacher's ALTS.
################################################################################
main_sim_fnc = function(num_sims, beta, X, p, sigma_1, sigma_2, q, include_bacher = TRUE) {
  results = list()
  i = 1
  while (i <= num_sims) {
    tryit = try(sim_fnc_1(beta, X, p, sigma_1, sigma_2, q = q, include_bacher = include_bacher), silent = FALSE)
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

################################################################################
## complete_sim_fnc(num_sims, beta, X, p, sigma_1, sigma_2, q = 1.35) 
##  For each `p` in the vector `p`, perform the simulation study for `num_sims` simulations.
##
## Arguments:
##  num_sims: Number of simulations to perform.
##  beta: Known regression coefficients.
##  X: Design matrix for the regression.
##  p: A vector of proportions of clean data to test.
##  sigma_1: Standard deviation of the clean data.
##  sigma_2: Standard deviation of the contaminated data.
##  q: The value of q to use in the new ALTS.
##
## Outputs:
##  The output is a list containing the following:
##    $Bacher: Results from Bacher's ALTS, further containing:
##      $p: Mean estimates for `p` from Bacher's ALTS at each true `p`.
##      $sigma: Mean estimates for `sigma_1` from Bacher's ALTS at each true `p`.
##      $rmape: Mean MAPE values from Bacher's ALTS at each true `p`.
##    $New: Results from the new ALTS.
##      $q: Value of `q` used from the new ALTS.
##      $p: Mean estimates for `p` from the new ALTS at each true `p`.
##      $sigma: Mean estimates for `sigma_1` from the new ALTS at each true `p`.
##      $rmape: Mean MAPE values from Bacher's ALTS at each true `p`.
################################################################################
complete_sim_fnc = function(num_sims, beta, X, p, sigma_1, sigma_2, q = 1.35) {
  bacher_p = matrix(data = NA, nrow = 1, ncol = length(p))
  bacher_sigma = matrix(data = NA, nrow = 1, ncol = length(p))
  bacher_rmape = matrix(data = NA, nrow = 1, ncol = length(p))
  new_q = matrix(data = NA, nrow = 1, ncol = length(p))
  new_p = matrix(data = NA, nrow = 1, ncol = length(p))
  new_sigma = matrix(data = NA, nrow = 1, ncol = length(p))
  new_rmape = matrix(data = NA, nrow = 1, ncol = length(p))
  for (j in 1:length(p)) {
    results = main_sim_fnc(num_sims, beta, X, p[j], sigma_1, sigma_2, q, include_bacher = TRUE)
    bacher_p[1, j] = results$Bacher$p
    bacher_sigma[1, j] = results$Bacher$sigma
    bacher_rmape[1, j] = results$Bacher$rmape
    new_q[1, j] = results$New$q
    new_p[1, j] = results$New$p
    new_sigma[1, j] = results$New$sigma
    new_rmape[1, j] = results$New$rmape
    print(j)
  }
  bacher_results = list(p = bacher_p, sigma = bacher_sigma, rmape = bacher_rmape)
  new_results = list(p = new_p, sigma = new_sigma, rmape = new_rmape, q = new_q)
  return(list(Bacher = bacher_results, New = new_results))
}

################################################################################
## complete_sim_fnc_vary_q(num_sims, beta, X, p, sigma_1, sigma_2, q) 
##  For each `p` in the vector `p`, and each `q` in the vector `q`, perform the 
##  simulation study for `num_sims` simulations.
##
## Arguments:
##  num_sims: Number of simulations to perform.
##  beta: Known regression coefficients.
##  X: Design matrix for the regression.
##  p: A vector of proportions of clean data to test.
##  sigma_1: Standard deviation of the clean data.
##  sigma_2: Standard deviation of the contaminated data.
##  q: A vector of values of q to use in the new ALTS.
##
## Outputs:
##  The output is a list containing the following:
##    $Bacher: Results from Bacher's ALTS, further containing:
##      $p: Mean estimates for `p` from Bacher's ALTS at each true `p`.
##      $sigma: Mean estimates for `sigma_1` from Bacher's ALTS at each true `p`.
##      $rmape: Mean MAPE values from Bacher's ALTS at each true `p`.
##    $New: Results from the new ALTS.
##      $p: Mean estimates for `p` from the new ALTS at each true `p`.
##      $sigma: Mean estimates for `sigma_1` from the new ALTS at each true `p`.
##      $rmape: Mean MAPE values from Bacher's ALTS at each true `p`.
################################################################################
complete_sim_fnc_vary_q = function(num_sims, beta, X, p, sigma_1, sigma_2, q) {
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
