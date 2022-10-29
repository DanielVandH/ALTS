################################################################################
## ALTS functions
################################################################################
## Contents (in order of presentation):
##  1: modified_mad
##          Compute the modified MAD estimator.
##
##  2. p_estimator(residual, sigma, q)
##          Compute the revised estimate of p = h/N.
##
##  3. weights_calculation(residual, p)
##          Compute the weights for LTS.
##
##  4. ALTS(model_formula, data, q = 1.0, maxIters = 20, tol = 1e-4, sigma_p = 0.25, residual_initial = rq(model_formula, data = data, method = "conquer")$residuals)
##          Fit model using ALTS.
##
##  5. ALTS_Bacher(model_formula, data, maxIters = 20, tol = 1e-4, residual_initial = rq(model_formula, data = data, method = "conquer")$residuals)
##          Fit model using Bacher's ALTS.
library(EnvStats)
library(L1pack)
library(dplyr)
library(caret)
library(formula.tools)
library(MLmetrics)
library(quantreg)
library(ggplot2)
library(pracma)

################################################################################
## modified_mad(residual, quantile)
##  Compute the modified MAD estimator.
##
## Arguments:
##  residual: Vector of residuals.
##  quantile: Estimate for p, the proportion of clean data.
##
## Outputs:
##  The modified MAD estimator, given as a standard deviation.
################################################################################
modified_mad = function(residual, quantile) {
  # Compute the numerator
  abs_e = abs(residual - median(residual))
  sort_e = sort(abs_e)
  h = floor(length(residual) * quantile)
  num_scale = sum(sort_e[1:h]^2)
  # Compute the denominator
  N = length(residual)
  p = (1:h)/(N+1)
  eta = qnormTrunc(p, mean = 0, sd = 1, min = 0, max = Inf)^2
  mean_eta = sum(eta)
  # Compute sigma 
  sigma = sqrt(num_scale/mean_eta)
  return(sigma)
}

################################################################################
## p_estimator(residual, sigma, q)
##  Compute the revised estimate of p = h/N.
##
## Arguments:
##  residual: Vector of residuals.
##  sigma: Estimate for the standard deviation of the clean data.
##  q: The value of `q` used in the ALTS algorithm.
##
## Outputs:
##  An estimate for `p`, the proportion of clean data.
################################################################################
p_estimator = function(residual, sigma, q) {
  # Compute E[Si^2]
  N = length(residual)
  p = (1:N)/(N+1)
  xi = qnormTrunc(p, mean = 0, sd = sigma, min = 0, max = Inf)
  Esi2 = cumsum(xi^2) 
  # Compute p 
  abs_e = abs(residual)
  sort_e = sort(abs_e)
  si2 = cumsum(sort_e^2)
  ratio = si2/Esi2 
  inequal = ratio <= q # ratio > q == FALSE  
  p = sum(inequal)/N
  return(p)
}

################################################################################
## weights_calculation(residual, p)
##  Compute the weights for LTS.
##
## Arguments:
##  residual: Vector of residuals.
##  p: An estimate for `p`, the proportion of clean data.
##
## Outputs:
##  The weights for LTS, given by `w_i = 1` if the absolute residual is less than
##    the `floor(Np)`th order statistic of the absolute residuals, and `0` otherwise.
################################################################################
weights_calculation = function(residual, p) {
  N = length(residual)
  h = floor(N*p)
  hc = N - h
  abs_e = abs(residual)
  sort_e_idx = sort(abs_e, index.return = TRUE)$ix
  weights = array()
  weights[sort_e_idx] = c(array(1, h), array(0, hc))
  return(weights)
}

################################################################################
## ALTS(model_formula, data, q = 1.0, maxIters = 20, tol = 1e-4, sigma_p = 0.25, 
##      residual_initial = rq(model_formula, data = data, method = "conquer")$residuals)
##  Fit model using ALTS.
##
## Arguments:
##  model_formula: The formula for the regression.
##  data: Data frame to use for the regression. Contains all the terms in model_formula.
##  q = 1.0: Value for the hyperparameter `q`.
##  maxIters = 20: The maximum allowable number of iterations.
##  tol = 1e-4: Threshold for the convergence, based on `abs(p[i] - p[i-1])/p[i] < tol`.
##  sigma_p = 0.25: Upper index to use in the modified MAD estimator.
##  residual_initial = rq(model_formula, data = data, method = "conquer")$residuals: Initial vector of residuals.
##
## Outputs:
##  The output is a list containing the following:
##    $p: The vector of p iterates.
##    $sigma: The vector of sigma iterates.
##    $beta: The final estimate of the regression coefficients.
##    $model: The final model.
##    $q: The q value used. This makes some of the simulation study parts a bit easier.
################################################################################
ALTS = function(model_formula, data, q = 1.0, maxIters = 20, tol = 1e-4, sigma_p = 0.25, 
                residual_initial = rq(model_formula, data = data, method = "conquer")$residuals) {
  sigma_initial = modified_mad(residual_initial, sigma_p)
  p_initial = 0.5
  data$weights = weights_calculation(residual_initial, p_initial)
  model_updated = glm(model_formula, weights = weights, data = data)
  # Iterate
  for (i in 2:maxIters) {
    updated_residuals = model_updated$residuals 
    sigma_initial[i] = modified_mad(updated_residuals, sigma_p)
    p_initial[i] = p_estimator(updated_residuals, sigma_initial[i], q)
    if (p_initial[i] < 0.5) {
      p_initial[i] = 0.5
      break
    }
    data$weights = weights_calculation(updated_residuals, p_initial[i])
    model_updated = glm(model_formula, weights = weights, data = data)
    if (abs(p_initial[i] - p_initial[i-1])/p_initial[i] < tol) {
      break
    }
  }
  # Return results 
  return(list(p = p_initial, sigma = sigma_initial, beta = model_updated$coefficients, model = model_updated, q = q))
}

################################################################################
## ALTS_Bacher(model_formula, data, maxIters = 20, tol = 1e-4, 
##      residual_initial = rq(model_formula, data = data, method = "conquer")$residuals)
##  Fit model using Bacher's ALTS.
##
## Arguments:
##  model_formula: The formula for the regression.
##  data: Data frame to use for the regression. Contains all the terms in model_formula.
##  maxIters = 20: The maximum allowable number of iterations.
##  tol = 1e-4: Threshold for the convergence, based on `abs(p[i] - p[i-1])/p[i] < tol`.
##  residual_initial = rq(model_formula, data = data, method = "conquer")$residuals: Initial vector of residuals.
##
## Outputs:
##  The output is a list containing the following:
##    $p: The vector of p iterates.
##    $sigma: The vector of sigma iterates.
##    $beta: The final estimate of the regression coefficients.
##    $model: The final model.
################################################################################
ALTS_Bacher = function(model_formula, data, maxIters = 20, tol = 1e-4, 
                       residual_initial = rq(model_formula, data = data, method = "conquer")$residuals) {
  N = nrow(data)
  MAR = median(abs(residual_initial - median(residual_initial)))
  sigma = MAR/qnormTrunc(0.5, mean = 0, sd = 1, min = 0, max = Inf)
  si2 = cumsum(sort((residual_initial)^2))/(1:N)
  p = sum((si2 < sigma^2))/N 
  data$weights = weights_calculation(residual_initial, p)
  model_updated = glm(model_formula, weights = weights, data = data)
  # Iterate
  sigma2_update = sigma^2 
  for (i in 2:maxIters) {
    h = floor(p[i-1]*N)
    sigma2_update[i] = mean(sort((model_updated$residuals)^2)[1:h])
    si2 = cumsum(sort((model_updated$residuals)^2))/(1:N)
    p[i] = sum((si2 < sigma2_update[i]))/N
    data$weights = weights_calculation(model_updated$residuals, p[i])
    model_updated = glm(model_formula, weights = weights, data = data)
    if (abs(p[i] - p[i-1])/p[i] < tol) {
      break
    }
  }
  # Return results 
  return(list(p = p, sigma = sqrt(sigma2_update), beta = model_updated$coefficients, model = model_updated))
}
