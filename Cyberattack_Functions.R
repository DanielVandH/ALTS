library(ggplot2)
library(magrittr)

## random_attack(mu, sigma, yt, p)
# Randomly attacks a proportion p of the data yt by scaling it by (1 + s/100), where s ~ Normal(mu, sigma^2).
random_attack = function(mu, sigma, yt, p) {
  prop = floor(p*length(yt))
  prop_idx = sample(1:length(yt), size = prop)
  s = rnorm(prop, mu, sigma)
  yt[prop_idx] = (1+s/100)*yt[prop_idx]
  return(yt)
}

## ramp_attack(t, yt, p, lambda, L)
# Randomly attacks a proportion p of the data yt using a ramp attack. 
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