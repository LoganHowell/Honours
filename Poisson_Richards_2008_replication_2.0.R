library(tidyverse)
library(Rfast)
library(epitools)
library(extraDistr)
library(foreach)
library(doParallel)
library(parallel)
library(stats4)
library(fAsianOptions)
library(maxLik)
library(matrixStats)
library(MASS)



tau = 10 # Number of minutes of observation. 
shade_stat = 0 # The status of the shading exposure, st. 0 = non-shaded, 1 = shaded. 
alpha = 0.5 # The arrival rate for one-flowered plants. 
alpha_shade = -0.05 # The constant component on the dummy variable associated with shaded plants. 
phi = 0.1 # The over-dispersion factor. 
r = 10 # The number of plants per 
num_reps = 600 # The number of model fits to be performed. 
Z = 10000 # The number of new points to use in calculation of cumulative IC (for each treatment). 
AIC_threshold = 6 # The threshold for delta_AIC-based selection. 


# Constraint sets for later optimisation:
# 2, 3, 2, 3

A_LL_M1 = matrix(c(1), 1, 1)
B_LL_M1 = matrix(c(0), 1, 1)

A_LL_M2 = matrix(c(1,0,0,1), 2, 2)
B_LL_M2 = matrix(c(0, .Machine$integer.max), 2, 1)

A_LL_M3 = matrix(c(1), 1, 1)
B_LL_M3 = matrix(c(0), 1, 1)

A_LL_M4 = matrix(c(1,1,0,1), 2, 2)
B_LL_M4 = matrix(c(0, 0), 2, 1)

A_LL_M5 = matrix(c(1,0,0,1), 2, 2)
B_LL_M5 = matrix(c(0, 0), 2, 1)

A_LL_M6 = matrix(c(1,0,0,0,1,0,0,0,1), 3,3)
B_LL_M6 = matrix(c(0, .Machine$integer.max, 0), 3, 1)

A_LL_M7 = matrix(c(1,0,0,1), 2, 2)
B_LL_M7 = matrix(c(0, 0), 2, 1)

A_LL_M8 = matrix(c(1,1,0,0,1,0,0,0,1), 3,3)
B_LL_M8 = matrix(c(0, 0, 0), 3, 1)


simplicity_key = matrix(
  data = c(
  c(0,1,1,1,1,1,1,1),
  c(0,0,0,1,1,1,1,1),
  c(1,1,0,1,1,1,1,1),
  c(0,1,0,0,1,1,1,1),
  c(0,1,0,1,0,1,1,1),
  c(0,0,0,0,0,0,0,1),
  c(0,1,0,1,1,1,0,1),
  c(0,0,0,0,0,1,0,0)
), nrow = 8, byrow = TRUE)



LL_M1 = function(pars, x)
{
  return(sum(x * log(pars[1] * tau) - pars[1] * tau - log(factorial(x))))
}

# @param: shade_stat is a vector of binary indicator variables determining whether
# the plant is in a shaded spot or not. 
LL_M2 = function(pars, x)
{
  return(sum(x * log((pars[1] + pars[2] * shade_stat) * tau) - (pars[1] + pars[2] * shade_stat) * tau - log(factorial(x))))
}

LL_M3 = function(pars, x)
{
  return(sum(x * log(pars[1] * x * tau) - pars[1] * x * tau - log(factorial(x))))  
}

LL_M4 = function(pars, x)
{
  return(sum(x * log((pars[1] + pars[2] * shade_stat) * x * tau) - 
               ((pars[1] + pars[2] * shade_stat) * x * tau) - log(factorial(x))))
}

# Check complete:
LL_M5 = function(pars, x)
{
  a = pars[1] / pars[2]
  b = 1 / pars[2]
  
  return(sum(lgamma(x + a) - lgamma(x + 1) - lgamma(a) + a * log((b / tau) / (1 + b/tau)) + x * log(1 / (1 + (b / tau)))))
}


LL_M6 = function(pars, x)
{
  a = (pars[1] + pars[2] * shade_stat) / pars[3]
  b = 1 / pars[3]
   
  return(sum(lgamma(x + a) - lgamma(x + 1) - lgamma(a) + a * log((b / tau) / (1 + b/tau)) + x * log(1 / (1 + (b / tau)))))  
}


LL_M7 = function(pars, x)
{
  a = (pars[1] * x) / pars[2]
  b = 1 / pars[2]

  return(sum(lgamma(x + a) - lgamma(x + 1) - lgamma(a) + a * log((b / tau) / (1 + b / tau)) + x * log(1 / (1 + (b / tau)))))
}


LL_M8 = function(pars, x)
{
  a = ((pars[1] + pars[2] * shade_stat) * x) / pars[3]
  b = 1 / pars[3]
  
  return(sum(lgamma(x + a) - lgamma(x + 1) - lgamma(a) + a * log((b / tau) / (1 + b / tau)) + x * log(1 / (1 + (b / tau)))))
}


flower_number = rep(1:4, each = r, times = 2) # Generate the selected flower groups.
shade_stat = rep(0:1, each = r * 4) # Generate the shade component. 

# We hold lambda constant as this reflects what is done within the paper: 
lambda_true = (alpha + alpha_shade * shade_stat) * flower_number ^ beta


process_rep = function(i)
{
  
  # Generate the i'th "truth" dataset: 
  truth_set = rep(0, times = length(flower_number))
  
  # Prevent the optimisation routines trying to incorporate lgamma(0) or log(0) and then failing: 
  while(0 %in% truth_set)
  {
    for(j in 1:length(flower_number))
    {
      b = 1 / phi
      truth_set[j] = stats::rnbinom(n = 1, size = lambda_true[j] / phi, prob = (b / tau) / (1 + b / tau))
    }
  }

  print(i)

  M1_fit = maxLik::maxLik(logLik = LL_M1, start = c(1), constraints = list(ineqA = A_LL_M1, ineqB = B_LL_M1), method = "BFGS", x = truth_set)
  M2_fit = maxLik::maxLik(logLik = LL_M2, start = c(1,1), constraints = list(ineqA = A_LL_M2, ineqB = B_LL_M2), x = truth_set)
  M3_fit = maxLik::maxLik(logLik = LL_M3, start = c(1), constraints = list(ineqA = A_LL_M3, ineqB = B_LL_M3), x = truth_set, method = "BFGS")
  M4_fit = maxLik::maxLik(logLik = LL_M4, start = c(1,1), constraints = list(ineqA = A_LL_M4, ineqB = B_LL_M4), x = truth_set)
  M5_fit = maxLik::maxLik(logLik = LL_M5, start = c(1,1), constraints = list(ineqA = A_LL_M5, ineqB = B_LL_M5), x = truth_set)
  M6_fit = maxLik::maxLik(logLik = LL_M6, start = c(1,1,1), constraints = list(ineqA = A_LL_M6, ineqB = B_LL_M6), x = truth_set)
  M7_fit = maxLik::maxLik(logLik = LL_M7, start = c(1,1), constraints = list(ineqA = A_LL_M7, ineqB = B_LL_M7), x = truth_set)
  M8_fit = maxLik::maxLik(logLik = LL_M8, start = c(1,1,1), constraints = list(ineqA = A_LL_M8, ineqB = B_LL_M8), x = truth_set)
  
  AIC_results = rep(0, times = 8)
  
  AIC_results[1] = -2 * M1_fit$maximum + 2 * length(M1_fit$estimate)
  AIC_results[2] = -2 * M2_fit$maximum + 2 * length(M2_fit$estimate)
  AIC_results[3] = -2 * M3_fit$maximum + 2 * length(M3_fit$estimate)
  AIC_results[4] = -2 * M4_fit$maximum + 2 * length(M4_fit$estimate)
  AIC_results[5] = -2 * M5_fit$maximum + 2 * length(M5_fit$estimate)
  AIC_results[6] = -2 * M6_fit$maximum + 2 * length(M6_fit$estimate)
  AIC_results[7] = -2 * M7_fit$maximum + 2 * length(M7_fit$estimate)
  AIC_results[8] = -2 * M8_fit$maximum + 2 * length(M8_fit$estimate)
  
  zth_results = rep(0, times = 8)
  
  for(z in 1:Z)
  {
    
    # Generate the i'th "truth" dataset: 
    x_valid = rep(0, times = length(flower_number))
    x_prob = rep(0, times = length(flower_number))
    
    # Prevent the optimisation routines trying to incorporate lgamma(0) or log(0) and then failing: 
 
    b =  1 / phi 
    
    x_valid = stats::rnbinom(n = length(flower_number), size = lambda_true / phi, prob = (b / tau) / (1 + b / tau))
    x_prob = stats::dnbinom(x = x_valid, size = lambda_true / phi, prob = (b / tau) / (1 + b / tau))
    
    c = sum((x_prob) * log(x_prob))
    
    # Calculate the KLD values themselves: 
    zth_results[1] = zth_results[1] + sum(log(x_prob)) - sum(log(dpois(x = x_valid, lambda = M1_fit$estimate[1] * tau)))
    
    
    zth_results[2] = zth_results[2] + sum(log(x_prob)) - sum(log(dpois(x = x_valid, lambda = (M2_fit$estimate[1] + M2_fit$estimate[2] * shade_stat) * tau)))
    
    
    zth_results[3] = zth_results[3] + sum(log(x_prob)) - sum(log(dpois(x = x_valid, lambda = M3_fit$estimate[1] * x_valid * tau)))
    
    
    zth_results[4] = zth_results[4] + sum(log(x_prob)) - sum(log(dpois(x = x_valid, lambda = (M4_fit$estimate[1] + M4_fit$estimate[2] * shade_stat) * x * tau)))
    
    
    size_5 = (M5_fit$estimate[1] * tau) / M5_fit$estimate[2]
    b_5 = 1 / M5_fit$estimate[2]
    prob_5 = (b_5 / tau) / (1 + b_5 / tau) 

    zth_results[5] = zth_results[5] + sum(log(x_prob)) - sum(log(stats::dnbinom(x = x_valid, size = size_5, prob = prob_5)))


    size_6 = ((M6_fit$estimate[1] + M6_fit$estimate[2] * shade_stat) * tau) / M6_fit$estimate[3]
    b_6 = 1 / M6_fit$estimate[3]
    prob_6 = (b_6 / tau) / (1 + b_6 / tau) 
    
    zth_results[6] = zth_results[6] + sum(log(x_prob)) - sum(log(stats::dnbinom(x = x_valid, size = size_6, prob = prob_6)))
    
    
    size_7 = (M7_fit$estimate[1] * x * tau) / M7_fit$estimate[2]
    b_7 = 1 / M7_fit$estimate[2]
    prob_7 = (b_7 / tau) / (1 + b_7 / tau) 
    
    zth_results[7] = zth_results[7] + sum(log(x_prob)) - sum(log(stats::dnbinom(x = x_valid, size = size_7, prob = prob_7)))
    
    
    size_8 = ((M8_fit$estimate[1] + M8_fit$estimate[2] * shade_stat) * x * tau) / M8_fit$estimate[3]
    b_8 = 1 / M8_fit$estimate[3]
    prob_8 = (b_8 / tau) / (1 + b_8 / tau) 
    
    zth_results[8] = zth_results[8] + sum(log(x_prob)) - sum(log(stats::dnbinom(x = x_valid, size = size_8, prob = prob_8)))
  }

  
  return_list = vector("list", length = 3) # Initialise a list to return from the function.
  
  # Save the return values in a set order:
  # @pos1: the mean KLD for the current fit.
  # @pos2: the AIC estimate for the current fit, calculated using the mean KLD values.
  # @pos3: the AIC values for the current fit, calculated using model max log likelihood.
  
  return_list[[1]] = zth_results / Z
  return_list[[2]] = 2 * (zth_results / Z - c)
  return_list[[3]] = AIC_results
  
  return(return_list)
}


# Set up the parallel environment:
num_cores = parallel::detectCores()
cluster = parallel::makeCluster(num_cores, type = "PSOCK")
registerDoSNOW(cluster)

pb = txtProgressBar(max = num_reps, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

parallel_results <- foreach(i=1:num_reps, 
                            .options.snow = opts) %dopar% process_rep(i)

close(pb)
stopCluster(cluster)

KLD_parallel_results = purrr::map(parallel_results, 1)
AIC_estimate_parallel_results = purrr::map(parallel_results, 2)
AIC_true_parallel_results = purrr::map(parallel_results, 3)


KLD_mat = matrix(unlist(KLD_parallel_results), ncol = 8, byrow = TRUE)
AIC_estimate_mat = matrix(unlist(AIC_estimate_parallel_results), ncol = 8, byrow = TRUE)
AIC_mat = matrix(unlist(AIC_true_parallel_results), ncol = 8, byrow = TRUE)



# Perform model selection based on the AIC values calculated above. 

# Calculate the AIC delta values. 
AIC_mat_delta = AIC_mat - rowMins(AIC_mat)
AIC_estimate_mat_delta = AIC_estimate_mat - rowMins(AIC_estimate_mat)

models_selected = list() # Create a list of models selected from each fit. 

for(i in 1:nrow(AIC_mat_delta))
{
  current_row = AIC_mat_delta[i, ] # Extract the i'th row. 
  
  # Determine which models are below the delta-AIC threshold. 
  models_selected[[i]] = which(current_row < AIC_threshold) 
  
  AIC_values = current_row[c(models_selected[[i]])]
  
  ordered_deltas = models_selected[[i]][order(AIC_values, decreasing = FALSE)]
  final_list = ordered_deltas
  
  for(j in 1:length(ordered_deltas))
  {
    removal_list = intersect(which(simplicity_key[ordered_deltas[j], ] == 1), ordered_deltas[j:length(ordered_deltas)])
    
    # If a more complicated model has a greater AIC value:
    if(length(removal_list) > 0)
    {
      final_list = final_list[!final_list %in% removal_list]
    }
  }
  
  # Overwrite the original list with the updated list (nested models removed)
  models_selected[[i]] = final_list
}

# Create an empty array to store the number of times each model is selected via the threshold selection rule. 
selection_times = rep(0, times = 8)

for(i in 1:length(models_selected))
{
  # Increment the selection counts for the models selected in the i'th fit:
  selection_times[models_selected[[i]]] = selection_times[models_selected[[i]]] + 1
}

# Calculate the proportion of times each model was selected. 
selection_prop = selection_times / num_reps

# Calculate mean KLD and the associated standard errors. 
EKLD = data.frame(
  EKLD = colMeans(KLD_mat),
  Model = c(1:8),
  SE_EKLD = matrixStats::colSds(KLD_mat) / sqrt(nrow(KLD_mat))
)

EKLD[order(EKLD$EKLD), ]