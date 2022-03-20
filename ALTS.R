library(EnvStats)
library(L1pack)
library(dplyr)
library(caret)

## Define the intermediate functions 

## modified_mad: Compute the modified MAD estimator.
modified_mad = function(residual, quantile) {
  # Compute the numerator
  abs_e = abs(residual - median(residual))
  sort_e = sort(abs_e)
  h = floor(length(residual) * quantile)
  num_scale = mean(sort_e[1:h])
  # Compute the denominator
  N = length(residual)
  p = (1:h)/(N+1)
  eta = qnormTrunc(p, mean = 0, sd = 1, min = 0, max = Inf)
  mean_eta = mean(eta)
  # Compute sigma 
  sigma = num_scale/mean_eta
  return(sigma)
}

## p_estimator: Compute the revised estimate of p = h/N.
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

## weights_calculation: Compute the weights for LTS. 
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

## ALTS: Fit model using ALTS 
ALTS = function(model_formula, data, q = 1.0, maxIters = 100, tol = 1e-4, use_lad = TRUE) {
  # Fit initial model and extract residuals 
  if (use_lad) {
    model_initial = lad(model_formula, data = data)
    residual_initial = model_initial$residuals
  } else {
    model_initial = rq(model_formula, data = data, method = "conquer")
    residual_initial = model_initial$residuals
  }
  # Update the model 
  sigma_initial = modified_mad(residual_initial, 0.5)
  p_initial = p_estimator(residual_initial, sigma_initial, q)
  if (p_initial < 0.5) {
    p_initial = 0.5
    return(list(p = p_initial, sigma = sigma_initial, beta = model_initial$coefficients, model = model_initial))
  }
  data$weights = weights_calculation(residual_initial, p_initial)
  model_updated = glm(model_formula, weights = weights, data = data)
  # Iterate
  for (i in 2:maxIters) {
    updated_residuals = model_updated$residuals 
    sigma_initial[i] = modified_mad(updated_residuals, p_initial[i-1])
    p_initial[i] = p_estimator(updated_residuals, sigma_initial[i], q)
    if (p_initial[i] < 0.5) {
      p_initial[i] = p_initial[i-1]
      break
    }
    data$weights = weights_calculation(updated_residuals, p_initial[i])
    model_updated = glm(model_formula, weights = weights, data = data)
    if (abs(p_initial[i] - p_initial[i-1])/p_initial[i] < tol) {
      break
    }
  }
  # Return results 
  return(list(p = p_initial, sigma = sigma_initial, beta = model_updated$coefficients, model = model_updated))
}

## ALTS_Bacher: Fit model using Bacher's ALTS
ALTS_Bacher = function(model_formula, data, maxIters = 100, tol = 1e-4, use_lad = TRUE) {
  N = nrow(data)
  # Fit initial model and extract residual
  if (use_lad) {
    model_initial = lad(model_formula, data = data)
    residual_initial = model_initial$residuals
  } else {
    model_initial = rq(model_formula, data = data, method = "conquer")
    residual_initial = model_initial$residuals
  }
  # Update the model
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

## Define method for ALTS in caret
lpALTS = list(type = "Regression",
              library = c("EnvStats", "L1pack", "dplyr"),
              loop = NULL, 
              prob = NULL,
              levels = NULL)

prm = data.frame(parameter = "q",
                 class = "numeric",
                 label = "Ratio")
lpALTS$parameters = prm 

altsGrid = function(x, y, len = NULL, search = "grid") {
  if (search == "grid") {
    out = data.frame(q = seq(1.0, 1.2, length.out = len))
  } else {
    out = data.frame(q = runif(len, min = 1.0, max = 1.2))
  }
  return(out)
}
lpALTS$grid = altsGrid 

altsFit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  model_formula = y ~ . 
  data = as.data.frame(x)
  data$y = y
  alts = ALTS(model_formula, data, q = param$q, ...)
  return(alts)
}
lpALTS$fit = altsFit

altsPred = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
  beta0 = modelFit$model$coefficients[1]
  beta1p = modelFit$model$coefficients[2:(length(modelFit$model$coefficients)-1)]
  mu = beta0 + newdata %*% beta1p 
  return(mu)
}
lpALTS$predict = altsPred 


altsSort = function(x) {
  x[order(x$q),]
}
lpALTS$sort = altsSort 


