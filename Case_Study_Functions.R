################################################################################
## Case study functions
################################################################################
## Contents (in order of presentation):
##  1: random_attack(mu, sigma, yt, p)
##          Randomly attacks a proportion `p` of the data `yt` by scaling it by `(1 + s/100)`, where `s ~ Normal(mu, sigma^2)`.
##
##  2. ramp_attack(t, yt, p, lambda, L)
##          Randomly attacks a proportion p of the data yt using a ramp attack. 
##
##  3. random_attack_ALTS(load_dat_training, load_dat_test, mu, sigma, p, q, model_formula)
##          Randomly attacks the data and then fits many models.
##
##  4. random_attack_all_simulations(num_sims, load_dat_training, load_dat_test, mu, sigma, p, q, model_formula)
##          Randomly attacks the data and then fits many models, repeated `num_sims` times.
##
##  5. random_attack_all_simulations_grid_vals(num_sims, load_dat_training, load_dat_test, mu, cv, p, q, model_formula)
##          Randomly attacks the data and then fits many models, repeated `num_sims` times, for a grid of values.
##
##  6. process_all_results(results, grid_vals)
##          Processes the results from the list of `results` from `random(ramp)_attack_all_simulations` that come from the corresponding parameters in each row of `grid_vals`.
##
##  7. ramp_attack_ALTS = function(load_dat_training, load_dat_test, lambda, L, p, q, model_formula)
##          Ramp attacks the data and then fits many models.
##
##  8. ramp_attack_all_simulations(num_sims, load_dat_training, load_dat_test, lambda, L, p, q, model_formula)
##          Ramp attacks the data and then fits many models, repeated `num_sims` times.
##
##  9. ramp_attack_all_simulations_grid_vals(num_sims, load_dat_training, load_dat_test, lambda, L, p, q, model_formula)
##          Ramp attacks the data and then fits many models, repeated `num_sims` times, for a grid of values.
##
##  10. ramp_attack_all_simulations_grid_vals_2(num_sims, load_dat_training, load_dat_test, gamma, L, p, q, model_formula)
##          Ramp attacks the data and then fits many models, repeated `num_sims` times, for a grid of values, based on the gamma formulation.
library(MASS)
library(rlmDataDriven)
library(ggplot2)

################################################################################
## random_attack(mu, sigma, yt, p)
##  Randomly attacks a proportion `p` of the data `yt` by scaling it by `(1 + s/100)`, 
##  where `s ~ Normal(mu, sigma^2)`.
##
## Arguments:
##  mu: The mean of `s`.
##  sigma: The standard deviation of `s`.
##  yt: The data to be attacked.
##  p: The proportion of data to be attacked.
##
## Outputs:
##  The attacked vector `yt`.
################################################################################
random_attack = function(mu, sigma, yt, p) {
  prop = floor(p*length(yt))
  prop_idx = sample(1:length(yt), size = prop)
  s = rnorm(prop, mu, sigma)
  yt[prop_idx] = (1+s/100)*yt[prop_idx]
  return(yt)
}

################################################################################
## ramp_attack(t, yt, p, lambda, L)
##  Randomly attacks a proportion p of the data yt using a ramp attack. 
##
## Arguments:
##  t: Trend vector.
##  yt: The data to be attacked.
##  p: The proportion of data to be attacked.
##  lambda: The scale parameter of the attack.
##  L: The length of each attack.
##
## Outputs:
##  The attacked vector `yt`.
################################################################################
ramp_attack = function(t, yt, p, lambda, L) {
  N = length(yt)
  num_groups = floor(N/L)
  num_groups_attack = floor(N*p/L)
  select_groups = sample(1:num_groups, size = num_groups_attack)
  time_intervals = cut_number(t, num_groups)
  attack_intervals = (1:length(levels(time_intervals)))[(levels(time_intervals) %in% levels(time_intervals)[select_groups])]
  split_data = split(data.frame(t, yt), time_intervals)
  for (i in attack_intervals) {
    data_attacking = split_data[[i]]
    ts = data_attacking$t[1]
    te = data_attacking$t %>% rev %>% .[1]
    first_idx = (ts < data_attacking$t) & (data_attacking$t <= (ts + te)/2)
    second_idx = ((ts+te)/2 < data_attacking$t) & (data_attacking$t < te)
    data_attacking$yt[first_idx] = (1 + lambda*(data_attacking$t[first_idx] - ts))*data_attacking$yt[first_idx]
    data_attacking$yt[second_idx] = (1 + lambda*(te - data_attacking$t[second_idx]))*data_attacking$yt[second_idx]
    split_data[[i]] = data_attacking
  }
  data = unsplit(split_data, time_intervals)$yt
}

################################################################################
## random_attack_ALTS = function(load_dat_training, load_dat_test, mu, sigma, p, q, model_formula)
##  Randomly attacks the data and then fits many models.
##
## Arguments:
##  load_dat_training: Training data to be used for the regression.
##  load_dat_test: Test data to be used for computing the MAPE.
##  mu: The mean of `s`.
##  sigma: The standard deviation of `s`.
##  p: The proportion of data to be attacked.
##  q: Value for the hyperparameter `q`.
##  model_formula: The formula for the regression.
##
## Outputs:
##  The result is a list with the following:
##    $New: Results from `ALTS`, the new ALTS method.
##    $Bacher: Results from `Bacher_ALTS`, Bacher's ALTS.
##    $Huber: Results for a Huber regression, fit using `rlm`.
##    $Bisquare: Results for a bisquare regression, fit using `rlm`.
##    $LeastSquares: Results for a least squares regression, fit using `lm`.
##    $Median: Results for a median regression, fit using `qr`.
##    $HuberDD: Results for a data-driven Huber regression, fit using `rlmDD`.
##    $BisquareDD: Results for a data-driven bisquare regression, fit using `rlmDD`.
################################################################################
random_attack_ALTS = function(load_dat_training, load_dat_test, mu, sigma, p, q, model_formula) {
  new_data = load_dat_training %>% mutate(Load = random_attack(mu, sigma, load_dat_training$Load, p))
  # New
  qr_results = rq(model_formula, data = new_data, method = "conquer")
  residual_initial = qr_results$residuals
  results = ALTS(model_formula, new_data, q, residual_initial = residual_initial)
  p_hat = results$p[length(results$p)]
  sigma_hat = results$sigma[length(results$sigma)]
  mape_val = MLmetrics::MAPE(predict(results$model, load_dat_test), load_dat_test$Load)
  new_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Bacher 
  results = ALTS_Bacher(model_formula, new_data, residual_initial = residual_initial)
  p_hat = results$p[length(results$p)]
  sigma_hat = results$sigma[length(results$sigma)]
  mape_val = MLmetrics::MAPE(predict(results$model, load_dat_test), load_dat_test$Load)
  bacher_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Huber benchmark
  huber_results = rlm(model_formula, new_data, psi = psi.huber)
  p_hat = NA
  sigma_hat = sigma(huber_results)
  mape_val = MLmetrics::MAPE(predict(huber_results, load_dat_test), load_dat_test$Load)
  huber_fixed_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Bisquare benchmark 
  bisquare_results = rlm(model_formula, new_data, psi = psi.bisquare)
  p_hat = NA
  sigma_hat = sigma(bisquare_results)
  mape_val = MLmetrics::MAPE(predict(bisquare_results, load_dat_test), load_dat_test$Load)
  bisquare_fixed_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Least squares benchmark 
  ls_results = lm(model_formula, new_data)
  p_hat = NA
  sigma_hat = sigma(ls_results)
  mape_val = MLmetrics::MAPE(predict(ls_results, load_dat_test), load_dat_test$Load)
  lsq_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Quantile regression
  p_hat = NA 
  sigma_hat = NA
  mape_val = MLmetrics::MAPE(predict(qr_results, load_dat_test), load_dat_test$Load)
  mqr_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Data-driven Huber
  y = new_data$Load
  x = model.matrix(model_formula, new_data)
  xx = x[,-1]
  huber_dd = rlmDD(y, xx, ls_results$coefficients, huber_results$coefficients, method = "Huber", plot = "N")
  p_hat = NA
  sigma_hat = NA
  mape_val = MLmetrics::MAPE(model.matrix(model_formula, load_dat_test)%*%huber_dd$esti$coefficients, load_dat_test$Load)
  HuberDD = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Data-driven bisquare
  y = new_data$Load
  x = model.matrix(model_formula, new_data)
  xx = x[,-1]
  bisquare_dd = rlmDD(y, xx, ls_results$coefficients, huber_results$coefficients, method = "Bisquare", plot = "N")
  p_hat = NA
  sigma_hat = NA
  mape_val = MLmetrics::MAPE(model.matrix(model_formula, load_dat_test)%*%bisquare_dd$esti$coefficients, load_dat_test$Load)
  BisquareDD = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Return
  return(list(New = new_results, Bacher = bacher_results, 
              Huber = huber_fixed_results, Bisquare = bisquare_fixed_results,
              LeastSquares = lsq_results, Median = mqr_results,
              HuberDD = HuberDD, BisquareDD = BisquareDD
  ))
}

################################################################################
## random_attack_all_simulations(num_sims, load_dat_training, load_dat_test, mu, sigma, p, q, model_formula)
##  Randomly attacks the data and then fits many models, repeated `num_sims` times.
##
## Arguments:
##  num_sims: Number of simulations to perform.
##  load_dat_training: Training data to be used for the regression.
##  load_dat_test: Test data to be used for computing the MAPE.
##  mu: The mean of `s`.
##  sigma: The standard deviation of `s`.
##  p: The proportion of data to be attacked.
##  q: Value for the hyperparameter `q`.
##  model_formula: The formula for the regression.
##
## Outputs:
##  The result is a list with the following:
##    $New: Results from `ALTS`, the new ALTS method.
##    $Bacher: Results from `Bacher_ALTS`, Bacher's ALTS.
##    $Huber: Results for a Huber regression, fit using `rlm`.
##    $Bisquare: Results for a bisquare regression, fit using `rlm`.
##    $LeastSquares: Results for a least squares regression, fit using `lm`.
##    $Median: Results for a median regression, fit using `qr`.
##    $HuberDD: Results for a data-driven Huber regression, fit using `rlmDD`.
##    $BisquareDD: Results for a data-driven bisquare regression, fit using `rlmDD`.
##  For each of these list elements we have the following sub-elements (for models where the estimate is not available, NA is used):
##    $p: An estimate for the proportion of clean data, `p`.
##    $sigma: An estimate for the standard deviation of the clean data, `sigma`.
##    $mape: The MAPE for the fitted model, computed on the test set.
##  These values are averaged over all simulations.
################################################################################
random_attack_all_simulations = function(num_sims, load_dat_training, load_dat_test, mu, sigma, p, q, model_formula) {
  new_p = c()
  new_sigma = c()
  new_mape = c()
  bacher_p = c()
  bacher_sigma = c()
  bacher_mape = c()
  huber_p = c()
  huber_sigma = c()
  huber_mape = c()
  bisquare_p = c()
  bisquare_sigma = c()
  bisquare_mape = c()
  ls_p = c()
  ls_sigma = c()
  ls_mape = c()
  qr_p = c()
  qr_sigma = c()
  qr_mape = c()
  huberdd_p = c()
  huberdd_sigma = c()
  huberdd_mape = c()
  bisquaredd_p = c()
  bisquaredd_sigma = c()
  bisquaredd_mape = c()
  for (i in 1:num_sims) {
    results = random_attack_ALTS(load_dat_training, load_dat_test, mu, sigma, p, q, model_formula)
    new_p = c(new_p, results$New$p)
    new_sigma = c(new_sigma, results$New$sigma)
    new_mape = c(new_mape, results$New$mape)
    bacher_p = c(bacher_p, results$Bacher$p)
    bacher_sigma = c(bacher_sigma, results$Bacher$sigma)
    bacher_mape = c(bacher_mape, results$Bacher$mape)
    huber_p = c(huber_p, results$Huber$p)
    huber_sigma = c(huber_sigma, results$Huber$sigma)
    huber_mape = c(huber_mape, results$Huber$mape)
    bisquare_p = c(bisquare_p, results$Bisquare$p)
    bisquare_sigma = c(bisquare_sigma, results$Bisquare$sigma)
    bisquare_mape = c(bisquare_mape, results$Bisquare$mape)
    ls_p = c(ls_p, results$LeastSquares$p)
    ls_sigma = c(ls_sigma, results$LeastSquares$sigma)
    ls_mape = c(ls_mape, results$LeastSquares$mape)
    qr_p = c(qr_p, results$Median$p)
    qr_sigma = c(qr_sigma, results$Median$sigma)
    qr_mape = c(qr_mape, results$Median$mape)
    huberdd_p = c(huberdd_p, results$HuberDD$p)
    huberdd_sigma = c(huberdd_sigma, results$HuberDD$sigma)
    huberdd_mape = c(huberdd_mape, results$HuberDD$mape)
    bisquaredd_p = c(bisquaredd_p, results$BisquareDD$p)
    bisquaredd_sigma = c(bisquaredd_sigma, results$BisquareDD$sigma)
    bisquaredd_mape = c(bisquaredd_mape, results$BisquareDD$mape)
  }
  new_p_mean = mean(new_p)
  new_sigma_mean = mean(new_sigma)
  new_mape_mean = mean(new_mape)
  new_results = list(p = new_p_mean, sigma = new_sigma_mean, mape = new_mape_mean)
  bacher_p_mean = mean(bacher_p)
  bacher_sigma_mean = mean(bacher_sigma)
  bacher_mape_mean = mean(bacher_mape)
  bacher_results = list(p = bacher_p_mean, sigma = bacher_sigma_mean, mape = bacher_mape_mean)
  huber_p_mean = mean(huber_p)
  huber_sigma_mean = mean(huber_sigma)
  huber_mape_mean = mean(huber_mape)
  huber_results = list(p = huber_p_mean, sigma = huber_sigma_mean, mape = huber_mape_mean)
  bisquare_p_mean = mean(bisquare_p)
  bisquare_sigma_mean = mean(bisquare_sigma)
  bisquare_mape_mean = mean(bisquare_mape)
  bisquare_results = list(p = bisquare_p_mean, sigma = bisquare_sigma_mean, mape = bisquare_mape_mean)
  ls_p_mean = mean(ls_p)
  ls_sigma_mean = mean(ls_sigma)
  ls_mape_mean = mean(ls_mape)
  ls_results = list(p = ls_p_mean, sigma = ls_sigma_mean, mape = ls_mape_mean)
  qr_p_mean = mean(qr_p)
  qr_sigma_mean = mean(qr_sigma)
  qr_mape_mean = mean(qr_mape)
  qr_results = list(p = qr_p_mean, sigma = qr_sigma_mean, mape = qr_mape_mean)
  huberdd_p_mean = mean(huberdd_p)
  huberdd_sigma_mean = mean(huberdd_sigma)
  huberdd_mape_mean = mean(huberdd_mape)
  huberdd_results = list(p = huberdd_p_mean, sigma = huberdd_sigma_mean, mape = huberdd_mape_mean)
  bisquaredd_p_mean = mean(bisquaredd_p)
  bisquaredd_sigma_mean = mean(bisquaredd_sigma)
  bisquaredd_mape_mean = mean(bisquaredd_mape)
  bisquaredd_results = list(p = bisquaredd_p_mean, sigma = bisquaredd_sigma_mean, mape = bisquaredd_mape_mean)
  return(list(New = new_results, Bacher = bacher_results, 
              Huber = huber_results, Bisquare = bisquare_results,
              LeastSquares = ls_results, Median = qr_results,
              HuberDD = huberdd_results, BisquareDD = bisquaredd_results
  ))
}

################################################################################
## random_attack_all_simulations_grid_vals(num_sims, load_dat_training, load_dat_test, mu, cv, p, q, model_formula)
##  Randomly attacks the data and then fits many models, repeated `num_sims` times, for a grid of values.
##
## Arguments:
##  num_sims: Number of simulations to perform.
##  load_dat_training: Training data to be used for the regression.
##  load_dat_test: Test data to be used for computing the MAPE.
##  mu: Vector of values for the mean of `s`.
##  cv: Vector for the coefficient of variations for `s`, defined by `cv = sigma/mu`.
##  p: Vector of proportions of data to be attacked.
##  q: Vector of values for the hyperparameter `q`.
##  model_formula: The formula for the regression.
##
## Outputs:
##  The result is a list with the following:
##    $results: A list of results from `random_attack_all_simulations`.
##    $grid_vals: Each row in `grid_vals` corresponds to the parameter values used 
##                to fit the corresponding element in `results`. This data frame should
##                have columns for `mu`, `cv`, `p`, and `q`.
################################################################################
random_attack_all_simulations_grid_vals = function(num_sims, load_dat_training, load_dat_test, mu, cv, p, q, model_formula) {
  grid_vals = expand.grid(mu = mu, cv = cv, p = p, q = q)
  results = list()
  for (i in 1:dim(grid_vals)[1]) {
    print(i)
    results[[i]] = random_attack_all_simulations(num_sims, load_dat_training, load_dat_test, grid_vals$mu[i], 
                                                 grid_vals$cv[i] * grid_vals$mu[i], grid_vals$p[i], grid_vals$q[i], model_formula) 
  }
  return(list(results = results, grid_vals = grid_vals))
}

################################################################################
## process_all_results(results, grid_vals)
##  Processes the results from the list of `results` from `random(ramp)_attack_all_simulations` that come from 
##  the corresponding parameters in each row of `grid_vals`.
##
## Arguments:
##  results: A list of results from `random(ramp)_attack_all_simulations` 
##  grid_vals: Each row in `grid_vals` corresponds to the parameter values used 
##             to fit the corresponding element in `results`.
##
## Outputs:
##  The result is a data frame with columns:
##    $phat: Estimates for the proportion of clean data, `p`.
##    $sigmahat: Estimates for the standard deviation of the clean data, `sigma`.
##    $mape: The MAPE for the fitted model, computed on the test set.
##    $Method: The method used for this row of results. See also `random_attack_all_simulations` for the methods used.
################################################################################
process_all_results = function(results, grid_vals) {
  new_p = c()
  new_sigma = c()
  new_mape = c()
  bacher_p = c()
  bacher_sigma = c()
  bacher_mape = c()
  huber_p = c()
  huber_sigma = c()
  huber_mape = c()
  bisquare_p = c()
  bisquare_sigma = c()
  bisquare_mape = c()
  ls_p = c()
  ls_sigma = c()
  ls_mape = c()
  qr_p = c()
  qr_sigma = c()
  qr_mape = c()
  huberdd_p = c()
  huberdd_sigma = c()
  huberdd_mape = c()
  bisquaredd_p = c()
  bisquaredd_sigma = c()
  bisquaredd_mape = c()
  for (i in 1:dim(grid_vals)[1]) {
    new_p = c(new_p, results[[i]]$New$p)
    new_sigma = c(new_sigma, results[[i]]$New$sigma)
    new_mape = c(new_mape, results[[i]]$New$mape)
    bacher_p = c(bacher_p, results[[i]]$Bacher$p)
    bacher_sigma = c(bacher_sigma, results[[i]]$Bacher$sigma)
    bacher_mape = c(bacher_mape, results[[i]]$Bacher$mape)
    huber_p = c(huber_p, results[[i]]$Huber$p)
    huber_sigma = c(huber_sigma, results[[i]]$Huber$sigma)
    huber_mape = c(huber_mape, results[[i]]$Huber$mape)
    bisquare_p = c(bisquare_p, results[[i]]$Bisquare$p)
    bisquare_sigma = c(bisquare_sigma, results[[i]]$Bisquare$sigma)
    bisquare_mape = c(bisquare_mape, results[[i]]$Bisquare$mape)
    ls_p = c(ls_p, results[[i]]$LeastSquares$p)
    ls_sigma = c(ls_sigma, results[[i]]$LeastSquares$sigma)
    ls_mape = c(ls_mape, results[[i]]$LeastSquares$mape)
    qr_p = c(qr_p, results[[i]]$Median$p)
    qr_sigma = c(qr_sigma, results[[i]]$Median$sigma)
    qr_mape = c(qr_mape, results[[i]]$Median$mape)
    huberdd_p = c(huberdd_p, results[[i]]$HuberDD$p)
    huberdd_sigma = c(huberdd_sigma, results[[i]]$HuberDD$sigma)
    huberdd_mape = c(huberdd_mape, results[[i]]$HuberDD$mape)
    bisquaredd_p = c(bisquaredd_p, results[[i]]$BisquareDD$p)
    bisquaredd_sigma = c(bisquaredd_sigma, results[[i]]$BisquareDD$sigma)
    bisquaredd_mape = c(bisquaredd_mape, results[[i]]$BisquareDD$mape)
  }
  new_grid_vals = grid_vals %>% mutate(phat = new_p, sigmahat = new_sigma, mape = new_mape, Method = "New")
  bacher_grid_vals = grid_vals %>% mutate(phat = bacher_p, sigmahat = bacher_sigma, mape = bacher_mape, Method = "Bacher")
  huber_grid_vals = grid_vals %>% mutate(phat = huber_p, sigmahat = huber_sigma, mape = huber_mape, Method = "Huber")
  bisquare_grid_vals = grid_vals %>% mutate(phat = bisquare_p, sigmahat = bisquare_sigma, mape = bisquare_mape, Method = "Bisquare")
  ls_grid_vals = grid_vals %>% mutate(phat = ls_p, sigmahat = ls_sigma, mape = ls_mape, Method = "LS")
  qr_grid_vals = grid_vals %>% mutate(phat = qr_p, sigmahat = qr_sigma, mape = qr_mape, Method = "Median")
  huberdd_grid_vals = grid_vals %>% mutate(phat = huberdd_p, sigmahat = huberdd_sigma, mape = huberdd_mape, Method = "Huber DD")
  bisquaredd_grid_vals = grid_vals %>% mutate(phat = bisquaredd_p, sigmahat = bisquaredd_sigma, mape = bisquaredd_mape, Method = "Bisquare DD")
  all_results = rbind(new_grid_vals, bacher_grid_vals, huber_grid_vals, bisquare_grid_vals, 
                      ls_grid_vals, qr_grid_vals, huberdd_grid_vals, bisquaredd_grid_vals)
  return(all_results)
}

################################################################################
## ramp_attack_ALTS = function(load_dat_training, load_dat_test, lambda, L, p, q, model_formula)
##  Ramp attacks the data and then fits many models.
##
## Arguments:
##  load_dat_training: Training data to be used for the regression.
##  load_dat_test: Test data to be used for computing the MAPE.
##  lambda: The scale parameter of the attack.
##  L: The length of each attack.
##  p: The proportion of data to be attacked.
##  q: Value for the hyperparameter `q`.
##  model_formula: The formula for the regression.
##
## Outputs:
##  The result is a list with the following:
##    $New: Results from `ALTS`, the new ALTS method.
##    $Bacher: Results from `Bacher_ALTS`, Bacher's ALTS.
##    $Huber: Results for a Huber regression, fit using `rlm`.
##    $Bisquare: Results for a bisquare regression, fit using `rlm`.
##    $LeastSquares: Results for a least squares regression, fit using `lm`.
##    $Median: Results for a median regression, fit using `qr`.
##    $HuberDD: Results for a data-driven Huber regression, fit using `rlmDD`.
##    $BisquareDD: Results for a data-driven bisquare regression, fit using `rlmDD`.
################################################################################
ramp_attack_ALTS = function(load_dat_training, load_dat_test, lambda, L, p, q, model_formula) {
  new_data = load_dat_training %>% mutate(Load = ramp_attack(load_dat_training$Trend, load_dat_training$Load, p, lambda, L))
  # New
  qr_results = rq(model_formula, data = new_data, method = "conquer")
  residual_initial = qr_results$residuals
  results = ALTS(model_formula, new_data, q, residual_initial = residual_initial)
  p_hat = results$p[length(results$p)]
  sigma_hat = results$sigma[length(results$sigma)]
  mape_val = MLmetrics::MAPE(predict(results$model, load_dat_test), load_dat_test$Load)
  new_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Bacher 
  results = ALTS_Bacher(model_formula, new_data, residual_initial = residual_initial)
  p_hat = results$p[length(results$p)]
  sigma_hat = results$sigma[length(results$sigma)]
  mape_val = MLmetrics::MAPE(predict(results$model, load_dat_test), load_dat_test$Load)
  bacher_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Huber benchmark
  huber_results = rlm(model_formula, new_data, psi = psi.huber)
  p_hat = NA
  sigma_hat = sigma(huber_results)
  mape_val = MLmetrics::MAPE(predict(huber_results, load_dat_test), load_dat_test$Load)
  huber_fixed_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Bisquare benchmark 
  bisquare_results = rlm(model_formula, new_data, psi = psi.bisquare)
  p_hat = NA
  sigma_hat = sigma(bisquare_results)
  mape_val = MLmetrics::MAPE(predict(bisquare_results, load_dat_test), load_dat_test$Load)
  bisquare_fixed_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Least squares benchmark 
  ls_results = lm(model_formula, new_data)
  p_hat = NA
  sigma_hat = sigma(ls_results)
  mape_val = MLmetrics::MAPE(predict(ls_results, load_dat_test), load_dat_test$Load)
  lsq_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Quantile regression
  p_hat = NA 
  sigma_hat = NA
  mape_val = MLmetrics::MAPE(predict(qr_results, load_dat_test), load_dat_test$Load)
  mqr_results = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Data-driven Huber
  y = new_data$Load
  x = model.matrix(model_formula, new_data)
  xx = x[,-1]
  huber_dd = rlmDD(y, xx, ls_results$coefficients, huber_results$coefficients, method = "Huber", plot = "N")
  p_hat = NA
  sigma_hat = NA
  mape_val = MLmetrics::MAPE(model.matrix(model_formula, load_dat_test)%*%huber_dd$esti$coefficients, load_dat_test$Load)
  HuberDD = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Data-driven bisquare
  y = new_data$Load
  x = model.matrix(model_formula, new_data)
  xx = x[,-1]
  bisquare_dd = rlmDD(y, xx, ls_results$coefficients, huber_results$coefficients, method = "Bisquare", plot = "N")
  p_hat = NA
  sigma_hat = NA
  mape_val = MLmetrics::MAPE(model.matrix(model_formula, load_dat_test)%*%bisquare_dd$esti$coefficients, load_dat_test$Load)
  BisquareDD = list(p = p_hat, sigma = sigma_hat, mape = mape_val)
  # Return
  return(list(New = new_results, Bacher = bacher_results, 
              Huber = huber_fixed_results, Bisquare = bisquare_fixed_results,
              LeastSquares = lsq_results, Median = mqr_results,
              HuberDD = HuberDD, BisquareDD = BisquareDD
  ))
}

################################################################################
## ramp_attack_all_simulations(num_sims, load_dat_training, load_dat_test, lambda, L, p, q, model_formula)
##  Ramp attacks the data and then fits many models, repeated `num_sims` times.
##
## Arguments:
##  num_sims: Number of simulations to perform.
##  load_dat_training: Training data to be used for the regression.
##  load_dat_test: Test data to be used for computing the MAPE.
##  lambda: The scale parameter of the attack.
##  L: The length of each attack.
##  p: The proportion of data to be attacked.
##  q: Value for the hyperparameter `q`.
##  model_formula: The formula for the regression.
##
## Outputs:
##  The result is a list with the following:
##    $New: Results from `ALTS`, the new ALTS method.
##    $Bacher: Results from `Bacher_ALTS`, Bacher's ALTS.
##    $Huber: Results for a Huber regression, fit using `rlm`.
##    $Bisquare: Results for a bisquare regression, fit using `rlm`.
##    $LeastSquares: Results for a least squares regression, fit using `lm`.
##    $Median: Results for a median regression, fit using `qr`.
##    $HuberDD: Results for a data-driven Huber regression, fit using `rlmDD`.
##    $BisquareDD: Results for a data-driven bisquare regression, fit using `rlmDD`.
##  For each of these list elements we have the following sub-elements (for models where the estimate is not available, NA is used):
##    $p: An estimate for the proportion of clean data, `p`.
##    $sigma: An estimate for the standard deviation of the clean data, `sigma`.
##    $mape: The MAPE for the fitted model, computed on the test set.
##  These values are averaged over all simulations.
################################################################################
ramp_attack_all_simulations = function(num_sims, load_dat_training, load_dat_test, lambda, L, p, q, model_formula) {
  new_p = c()
  new_sigma = c()
  new_mape = c()
  bacher_p = c()
  bacher_sigma = c()
  bacher_mape = c()
  huber_p = c()
  huber_sigma = c()
  huber_mape = c()
  bisquare_p = c()
  bisquare_sigma = c()
  bisquare_mape = c()
  ls_p = c()
  ls_sigma = c()
  ls_mape = c()
  qr_p = c()
  qr_sigma = c()
  qr_mape = c()
  huberdd_p = c()
  huberdd_sigma = c()
  huberdd_mape = c()
  bisquaredd_p = c()
  bisquaredd_sigma = c()
  bisquaredd_mape = c()
  for (i in 1:num_sims) {
    results = ramp_attack_ALTS(load_dat_training, load_dat_test, lambda, L, p, q, model_formula)
    new_p = c(new_p, results$New$p)
    new_sigma = c(new_sigma, results$New$sigma)
    new_mape = c(new_mape, results$New$mape)
    bacher_p = c(bacher_p, results$Bacher$p)
    bacher_sigma = c(bacher_sigma, results$Bacher$sigma)
    bacher_mape = c(bacher_mape, results$Bacher$mape)
    huber_p = c(huber_p, results$Huber$p)
    huber_sigma = c(huber_sigma, results$Huber$sigma)
    huber_mape = c(huber_mape, results$Huber$mape)
    bisquare_p = c(bisquare_p, results$Bisquare$p)
    bisquare_sigma = c(bisquare_sigma, results$Bisquare$sigma)
    bisquare_mape = c(bisquare_mape, results$Bisquare$mape)
    ls_p = c(ls_p, results$LeastSquares$p)
    ls_sigma = c(ls_sigma, results$LeastSquares$sigma)
    ls_mape = c(ls_mape, results$LeastSquares$mape)
    qr_p = c(qr_p, results$Median$p)
    qr_sigma = c(qr_sigma, results$Median$sigma)
    qr_mape = c(qr_mape, results$Median$mape)
    huberdd_p = c(huberdd_p, results$HuberDD$p)
    huberdd_sigma = c(huberdd_sigma, results$HuberDD$sigma)
    huberdd_mape = c(huberdd_mape, results$HuberDD$mape)
    bisquaredd_p = c(bisquaredd_p, results$BisquareDD$p)
    bisquaredd_sigma = c(bisquaredd_sigma, results$BisquareDD$sigma)
    bisquaredd_mape = c(bisquaredd_mape, results$BisquareDD$mape)
  }
  new_p_mean = mean(new_p)
  new_sigma_mean = mean(new_sigma)
  new_mape_mean = mean(new_mape)
  new_results = list(p = new_p_mean, sigma = new_sigma_mean, mape = new_mape_mean)
  bacher_p_mean = mean(bacher_p)
  bacher_sigma_mean = mean(bacher_sigma)
  bacher_mape_mean = mean(bacher_mape)
  bacher_results = list(p = bacher_p_mean, sigma = bacher_sigma_mean, mape = bacher_mape_mean)
  huber_p_mean = mean(huber_p)
  huber_sigma_mean = mean(huber_sigma)
  huber_mape_mean = mean(huber_mape)
  huber_results = list(p = huber_p_mean, sigma = huber_sigma_mean, mape = huber_mape_mean)
  bisquare_p_mean = mean(bisquare_p)
  bisquare_sigma_mean = mean(bisquare_sigma)
  bisquare_mape_mean = mean(bisquare_mape)
  bisquare_results = list(p = bisquare_p_mean, sigma = bisquare_sigma_mean, mape = bisquare_mape_mean)
  ls_p_mean = mean(ls_p)
  ls_sigma_mean = mean(ls_sigma)
  ls_mape_mean = mean(ls_mape)
  ls_results = list(p = ls_p_mean, sigma = ls_sigma_mean, mape = ls_mape_mean)
  qr_p_mean = mean(qr_p)
  qr_sigma_mean = mean(qr_sigma)
  qr_mape_mean = mean(qr_mape)
  qr_results = list(p = qr_p_mean, sigma = qr_sigma_mean, mape = qr_mape_mean)
  huberdd_p_mean = mean(huberdd_p)
  huberdd_sigma_mean = mean(huberdd_sigma)
  huberdd_mape_mean = mean(huberdd_mape)
  huberdd_results = list(p = huberdd_p_mean, sigma = huberdd_sigma_mean, mape = huberdd_mape_mean)
  bisquaredd_p_mean = mean(bisquaredd_p)
  bisquaredd_sigma_mean = mean(bisquaredd_sigma)
  bisquaredd_mape_mean = mean(bisquaredd_mape)
  bisquaredd_results = list(p = bisquaredd_p_mean, sigma = bisquaredd_sigma_mean, mape = bisquaredd_mape_mean)
  return(list(New = new_results, Bacher = bacher_results, 
              Huber = huber_results, Bisquare = bisquare_results,
              LeastSquares = ls_results, Median = qr_results,
              HuberDD = huberdd_results, BisquareDD = bisquaredd_results
  ))
}

################################################################################
## ramp_attack_all_simulations_grid_vals(num_sims, load_dat_training, load_dat_test, lambda, L, p, q, model_formula)
##  Ramp attacks the data and then fits many models, repeated `num_sims` times, for a grid of values.
##
## Arguments:
##  num_sims: Number of simulations to perform.
##  load_dat_training: Training data to be used for the regression.
##  load_dat_test: Test data to be used for computing the MAPE.
##  lambda: Vector of scale parameters.
##  L: Vector of lengths for each attack.
##  p: Vector of proportions of data to be attacked.
##  q: Vector of values for the hyperparameter `q`.
##  model_formula: The formula for the regression.
##
## Outputs:
##  The result is a list with the following:
##    $results: A list of results from `ramp_attack_all_simulations`.
##    $grid_vals: Each row in `grid_vals` corresponds to the parameter values used 
##                to fit the corresponding element in `results`. This data frame should
##                have columns for `lambda`, `L`, `p`, and `q`.
################################################################################
ramp_attack_all_simulations_grid_vals = function(num_sims, load_dat_training, load_dat_test, lambda, L, p, q, model_formula) {
  grid_vals = expand.grid(L = L, lambda = lambda, p = p, q = q)
  results = list()
  for (i in 1:dim(grid_vals)[1]) {
    print(i)
    results[[i]] = ramp_attack_all_simulations(num_sims, load_dat_training, load_dat_test, 
                                               grid_vals$lambda[i], grid_vals$L[i], grid_vals$p[i], grid_vals$q[i], model_formula) 
  }
  return(list(results = results, grid_vals = grid_vals))
}

################################################################################
## ramp_attack_all_simulations_grid_vals_2(num_sims, load_dat_training, load_dat_test, gamma, L, p, q, model_formula)
##  Ramp attacks the data and then fits many models, repeated `num_sims` times, for a grid of values, based on the gamma formulation.
##
## Arguments:
##  num_sims: Number of simulations to perform.
##  load_dat_training: Training data to be used for the regression.
##  load_dat_test: Test data to be used for computing the MAPE.
##  gamma: Vector of gamma parameters, related to the scale parameters by `lambda = 1 + gamma*L/2`.
##  L: Vector of lengths for each attack.
##  p: Vector of proportions of data to be attacked.
##  q: Vector of values for the hyperparameter `q`.
##  model_formula: The formula for the regression.
##
## Outputs:
##  The result is a list with the following:
##    $results: A list of results from `ramp_attack_all_simulations`.
##    $grid_vals: Each row in `grid_vals` corresponds to the parameter values used 
##                to fit the corresponding element in `results`. This data frame should
##                have columns for `gamma`, `L`, `p`, and `q`.
################################################################################
ramp_attack_all_simulations_grid_vals_2 = function(num_sims, load_dat_training, load_dat_test, gamma, L, p, q, model_formula) {
  grid_vals = expand.grid(L = L, gamma = gamma, p = p, q = q)
  results = list()
  for (i in 1:dim(grid_vals)[1]) {
    print(i)
    results[[i]] = ramp_attack_all_simulations(num_sims, load_dat_training, load_dat_test, 2*(grid_vals$gamma[i] - 1)/grid_vals$L[i], grid_vals$L[i], grid_vals$p[i], grid_vals$q[i], model_formula) 
  }
  return(list(results = results, grid_vals = grid_vals))
}

################################################################################
##  widen_results(data)
##    Widens the results from an experiment.
##
##  Arguments:
##    data: The data frame of results to widen.
##
##  Outputs:
##    The output is the widened `data`.
################################################################################
widen_results_and_latexify = function(data) {
  data %>%
    mutate(mape = 100*mape) %>%
    ungroup() %>% # Needed for selecting q 
    select(-phat, -q, -sigmahat) %>%
    pivot_wider(names_from = "Method", values_from = "mape") %>%
    .[, c(3, 1, 2, 5, 4, 7, 11, 6, 10, 9, 8)] %>% 
    mutate(p = 1 - p) %>% 
    mutate(across(Jiao:LS, function(x) sprintf("%.2f", round(x, digits = 2)))) %>% 
    kable(booktabs = TRUE)
}
