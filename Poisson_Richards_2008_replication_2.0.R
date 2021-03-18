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
library(ramify)
library(purrr)
library(doSNOW)
library(latex2exp)
library(ggrepel)
library(data.table)
library(plotly)
library(cowplot)
library(xlsx)

tau = 10 # Number of minutes of observation. 
shade_stat = 0 # The status of the shading exposure, st. 0 = non-shaded, 1 = shaded. 
alpha = 0.5 # The arrival rate for one-flowered plants. 
alpha_shade = -0.05 # The constant component on the dummy variable associated with shaded plants. 
beta = 0.8
phi = 0.1 # The over-dispersion factor. 
r = 5 # The number of plants per replicate. 
num_reps = 1000 # The number of model fits to be performed. 
Z = 10000 # The number of new points to use in calculation of cumulative IC (for each treatment). 
AIC_threshold = 6 # The threshold for delta_AIC-based selection. 
# v_tilda = 2 # Set the variance inflation factor manually. 


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



LL_M1 = function(pars, x_arrivals)
{
  
  lambda = pars[1]
  
  return(sum(-lambda * tau + x_arrivals * log(lambda * tau) - log(factorial(x_arrivals))))
}

# @param: shade_stat is a vector of binary indicator variables determining whether
# the plant is in a shaded spot or not. 
LL_M2 = function(pars, x_arrivals, shade_stat_vec)
{
  lambda = pars[1] + pars[2] * shade_stat_vec
  
  return(sum(-lambda * tau + x_arrivals * log(lambda * tau) - log(factorial(x_arrivals))))
}

LL_M3 = function(pars, x_arrivals, y_flower_number)
{
  
  lambda = pars[1] * y_flower_number
  
  return(sum(-lambda * tau + x_arrivals * log(lambda * tau) - log(factorial(x_arrivals))))  
}

LL_M4 = function(pars, x_arrivals, y_flower_number, shade_stat_vec)
{
  lambda = (pars[1] + pars[2] * shade_stat_vec) * y_flower_number
  
  return(sum(-lambda * tau + x_arrivals * log(lambda * tau) - log(factorial(x_arrivals))))
}


LL_M5 = function(pars, x_arrivals)
{
  lambda = pars[1]
  
  a = lambda / pars[2]
  b = 1 / pars[2]
  
  return(sum(lgamma(x_arrivals + a) - lgamma(x_arrivals + 1) - lgamma(a) + a * log((b / tau) / (1 + b / tau)) + x_arrivals * log(1 / (1 + (b / tau)))))
}


LL_M6 = function(pars, x_arrivals, shade_stat_vec)
{
  lambda = pars[1] + pars[2] * shade_stat_vec
  
  a = lambda / pars[3]
  b = 1 / pars[3]
   
  return(sum(lgamma(x_arrivals + a) - lgamma(x_arrivals + 1) - lgamma(a) + a * log((b / tau) / (1 + b / tau)) + x_arrivals * log(1 / (1 + (b / tau)))))
}


LL_M7 = function(pars, x_arrivals, y_flower_number)
{
  lambda = pars[1] * y_flower_number
  
  a = lambda / pars[2]
  b = 1 / pars[2]

  return(sum(lgamma(x_arrivals + a) - lgamma(x_arrivals + 1) - lgamma(a) + a * log((b / tau) / (1 + b / tau)) + x_arrivals * log(1 / (1 + (b / tau)))))
}


LL_M8 = function(pars, x_arrivals, y_flower_number, shade_stat_vec)
{
  lambda = (pars[1] + pars[2] * shade_stat_vec) * y_flower_number
  
  a = lambda / pars[3]
  b = 1 / pars[3]
  
  return(sum(lgamma(x_arrivals + a) - lgamma(x_arrivals + 1) - lgamma(a) + a * log((b / tau) / (1 + b / tau)) + x_arrivals * log(1 / (1 + (b / tau)))))
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
  #while(0 %in% truth_set)
  #{
    for(j in 1:length(flower_number))
    {
      b = 1 / phi
      truth_set[j] = stats::rnbinom(n = 1, size = lambda_true[j] / phi, prob = (b / tau) / (1 + b / tau))
    }
  #}

  M1_fit = maxLik::maxLik(logLik = LL_M1, start = c(1), constraints = list(ineqA = A_LL_M1, ineqB = B_LL_M1), method = "BFGS", x_arrivals = truth_set)
  M2_fit = maxLik::maxLik(logLik = LL_M2, start = c(1,1), constraints = list(ineqA = A_LL_M2, ineqB = B_LL_M2), x_arrivals = truth_set, shade_stat_vec = shade_stat)
  M3_fit = maxLik::maxLik(logLik = LL_M3, start = c(1), constraints = list(ineqA = A_LL_M3, ineqB = B_LL_M3), method = "BFGS", x_arrivals = truth_set, y_flower_number = flower_number)
  M4_fit = maxLik::maxLik(logLik = LL_M4, start = c(1,1), constraints = list(ineqA = A_LL_M4, ineqB = B_LL_M4), x_arrivals = truth_set, y_flower_number = flower_number, shade_stat_vec = shade_stat)
  M5_fit = maxLik::maxLik(logLik = LL_M5, start = c(1,1), constraints = list(ineqA = A_LL_M5, ineqB = B_LL_M5), x_arrivals = truth_set)
  M6_fit = maxLik::maxLik(logLik = LL_M6, start = c(1,1,1), constraints = list(ineqA = A_LL_M6, ineqB = B_LL_M6), x_arrivals = truth_set, shade_stat_vec = shade_stat)
  M7_fit = maxLik::maxLik(logLik = LL_M7, start = c(1,1), constraints = list(ineqA = A_LL_M7, ineqB = B_LL_M7), x_arrivals = truth_set, y_flower_number = flower_number)
  M8_fit = maxLik::maxLik(logLik = LL_M8, start = c(1,1,1), constraints = list(ineqA = A_LL_M8, ineqB = B_LL_M8), x_arrivals = truth_set, y_flower_number = flower_number, shade_stat_vec = shade_stat)
  
  AIC_results = rep(0, times = 8)
  
  AIC_results[1] = -2 * M1_fit$maximum + 2 * length(M1_fit$estimate)
  AIC_results[2] = -2 * M2_fit$maximum + 2 * length(M2_fit$estimate)
  AIC_results[3] = -2 * M3_fit$maximum + 2 * length(M3_fit$estimate)
  AIC_results[4] = -2 * M4_fit$maximum + 2 * length(M4_fit$estimate)
  AIC_results[5] = -2 * M5_fit$maximum + 2 * length(M5_fit$estimate)
  AIC_results[6] = -2 * M6_fit$maximum + 2 * length(M6_fit$estimate)
  AIC_results[7] = -2 * M7_fit$maximum + 2 * length(M7_fit$estimate)
  AIC_results[8] = -2 * M8_fit$maximum + 2 * length(M8_fit$estimate)
  
  
  # Perform the required QAIC calculations: 
  
  QAIC_results = rep(0, times = 8)
  
  saturated_log_likelihood = exp(-truth_set) * (truth_set ^ truth_set) / factorial(truth_set)
  
  saturated_log_likelihood = sum(log(saturated_log_likelihood))
  
  df = length(shade_stat) - 2
  v_tilda = (2 / df) * (saturated_log_likelihood - M4_fit$maximum)
  
  QAIC_results[1] = - (2 / v_tilda) * M1_fit$maximum + 2 * length(M1_fit$estimate)
  QAIC_results[2] = - (2 / v_tilda) * M2_fit$maximum + 2 * length(M2_fit$estimate)
  QAIC_results[3] = - (2 / v_tilda) * M3_fit$maximum + 2 * length(M3_fit$estimate)
  QAIC_results[4] = - (2 / v_tilda) * M4_fit$maximum + 2 * length(M4_fit$estimate)
  QAIC_results[5] = - (2 / v_tilda) * M5_fit$maximum + 2 * length(M5_fit$estimate)
  QAIC_results[6] = - (2 / v_tilda) * M6_fit$maximum + 2 * length(M6_fit$estimate)
  QAIC_results[7] = - (2 / v_tilda) * M7_fit$maximum + 2 * length(M7_fit$estimate)
  QAIC_results[8] = - (2 / v_tilda) * M8_fit$maximum + 2 * length(M8_fit$estimate)
  
  
  zth_results = rep(0, times = 8)
  
  c = 0 # Initialise a value to estimate the contstant, c. 
  
  for(z in 1:Z)
  {
    
    # Generate the i'th "truth" dataset: 
    x_valid = rep(0, times = length(flower_number))
    x_prob = rep(0, times = length(flower_number))
  
    b =  1 / phi 
    
    x_valid = stats::rnbinom(n = length(flower_number), size = lambda_true / phi, prob = (b / tau) / (1 + b / tau))
    x_prob = stats::dnbinom(x = x_valid, size = lambda_true / phi, prob = (b / tau) / (1 + b / tau))
    
    c = c + sum(log(x_prob))
    
    # Calculate the KLD values themselves: 
    
    M1_lambda = M1_fit$estimate[1]
    zth_results[1] = zth_results[1] + sum(log(x_prob)) - sum(log(dpois(x = x_valid, lambda = M1_lambda * tau)))
    
    
    M2_lambda = M2_fit$estimate[1] + M2_fit$estimate[2] * shade_stat
    zth_results[2] = zth_results[2] + sum(log(x_prob)) - sum(log(dpois(x = x_valid, lambda = M2_lambda * tau)))
    
    
    M3_lambda = M3_fit$estimate[1] * flower_number
    zth_results[3] = zth_results[3] + sum(log(x_prob)) - sum(log(dpois(x = x_valid, lambda = M3_lambda * tau)))
    
    
    M4_lambda = (M4_fit$estimate[1] + M4_fit$estimate[2] * shade_stat) * flower_number
    zth_results[4] = zth_results[4] + sum(log(x_prob)) - sum(log(dpois(x = x_valid, lambda = M4_lambda * tau)))
    
    
    size_5 = (M5_fit$estimate[1]) / M5_fit$estimate[2]
    b_5 = 1 / M5_fit$estimate[2]
    prob_5 = (b_5 / tau) / (1 + b_5 / tau) 

    zth_results[5] = zth_results[5] + sum(log(x_prob)) - sum(log(stats::dnbinom(x = x_valid, size = size_5, prob = prob_5)))


    size_6 = ((M6_fit$estimate[1] + M6_fit$estimate[2] * shade_stat)) / M6_fit$estimate[3]
    b_6 = 1 / M6_fit$estimate[3]
    prob_6 = (b_6 / tau) / (1 + b_6 / tau) 
    
    zth_results[6] = zth_results[6] + sum(log(x_prob)) - sum(log(stats::dnbinom(x = x_valid, size = size_6, prob = prob_6)))
    
    
    size_7 = (M7_fit$estimate[1] * flower_number) / M7_fit$estimate[2]
    b_7 = 1 / M7_fit$estimate[2]
    prob_7 = (b_7 / tau) / (1 + b_7 / tau) 
    
    zth_results[7] = zth_results[7] + sum(log(x_prob)) - sum(log(stats::dnbinom(x = x_valid, size = size_7, prob = prob_7)))
    
    
    size_8 = ((M8_fit$estimate[1] + M8_fit$estimate[2] * shade_stat) * flower_number) / M8_fit$estimate[3]
    b_8 = 1 / M8_fit$estimate[3]
    prob_8 = (b_8 / tau) / (1 + b_8 / tau) 
    
    zth_results[8] = zth_results[8] + sum(log(x_prob)) - sum(log(stats::dnbinom(x = x_valid, size = size_8, prob = prob_8)))
  }

  
  incorrect_inference = rep(0, length = 8)
  
  incorrect_inference[[2]] = as.integer(M2_fit$estimate[[2]] > 0)
  incorrect_inference[[4]] = as.integer(M4_fit$estimate[[2]] > 0)
  incorrect_inference[[6]] = as.integer(M6_fit$estimate[[2]] > 0)
  incorrect_inference[[8]] = as.integer(M8_fit$estimate[[2]] > 0)
  
  return_list = vector("list", length = 5) # Initialise a list to return from the function.
  
  # Save the return values in a set order:
  # @pos1: the mean KLD for the current fit.
  # @pos2: the AIC estimate for the current fit, calculated using the mean KLD values.
  # @pos3: the AIC values for the current fit, calculated using model max log likelihood.
  
  return_list[[1]] = zth_results / Z
  return_list[[2]] = 2 * (zth_results / Z - c / Z)
  return_list[[3]] = AIC_results
  return_list[[4]] = QAIC_results
  return_list[[5]] = incorrect_inference
  
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

KLD_parallel_results = purrr::map(parallel_results, 1)
AIC_estimate_parallel_results = purrr::map(parallel_results, 2)
AIC_true_parallel_results = purrr::map(parallel_results, 3)
QAIC_true_parallel_results = purrr::map(parallel_results, 4)
incorrect_inference_results = purrr::map(parallel_results, 5)

KLD_mat = matrix(unlist(KLD_parallel_results), ncol = 8, byrow = TRUE)
AIC_estimate_mat = matrix(unlist(AIC_estimate_parallel_results), ncol = 8, byrow = TRUE)
AIC_mat = matrix(unlist(AIC_true_parallel_results), ncol = 8, byrow = TRUE)
QAIC_mat = matrix(unlist(QAIC_true_parallel_results), ncol = 8, byrow = TRUE)
inference_matrix =  matrix(unlist(incorrect_inference_results), ncol = 8, byrow = TRUE)


xlsx::write.xlsx(AIC_mat, file = "poisson_scenario_f.xlsx", sheetName = "AIC", row.names = FALSE, col.names = FALSE)
xlsx::write.xlsx(KLD_mat, file = "poisson_scenario_f.xlsx", sheetName = "KLD_estimates", row.names = FALSE, col.names = FALSE, append = TRUE)
xlsx::write.xlsx(AIC_estimate_mat, file = "poisson_scenario_f.xlsx", sheetName = "AIC_KLD_estimates", row.names = FALSE, col.names = FALSE, append = TRUE)
xlsx::write.xlsx(QAIC_mat, file = "poisson_scenario_f.xlsx", sheetName = "QAIC", row.names = FALSE, col.names = FALSE, append = TRUE)
xlsx::write.xlsx(incorrect_inference_results, file = "poisson_scenario_f.xlsx", sheetName = "incorrect_inference", row.names = FALSE, col.names = FALSE, append = TRUE)






QAIC_mat = QAIC_mat[, 1:4]







colmeans(AIC_mat)
colmeans(AIC_estimate_mat)
colmeans(QAIC_mat)

# Perform model selection based on the AIC values calculated above. 

# Calculate the AIC delta values. 
AIC_mat_delta = AIC_mat - rowMins(AIC_mat)
AIC_estimate_mat_delta = AIC_estimate_mat - rowMins(AIC_estimate_mat)
QAIC_mat_delta = QAIC_mat - rowMins(QAIC_mat)




prob_model_selected = function(AIC_threshold, QAIC_threshold, AIC_mat_delta, QAIC_mat_delta)
{
  # Create a list of models selected from each fit.
  nested_AIC_models_selected = list() 
  delta_AIC_models_selected = list()
  delta_QAIC_models_selected = list()
  nested_QAIC_models_selected = list()
  
  for(i in 1:nrow(AIC_mat_delta))
  {
    # Extract the i'th rows. 
    current_AIC_row = AIC_mat_delta[i, ] 
    current_QAIC_row = QAIC_mat_delta[i, ]
    
    # Determine which models are below the delta-AIC threshold. 
    nested_AIC_models_selected[[i]] = which(current_AIC_row < AIC_threshold) 
    delta_AIC_models_selected[[i]] = which(current_AIC_row < AIC_threshold)
    delta_QAIC_models_selected[[i]] = which(current_QAIC_row < QAIC_threshold)
    nested_QAIC_models_selected[[i]] = which(current_QAIC_row < QAIC_threshold)
    
    AIC_values = current_AIC_row[c(nested_AIC_models_selected[[i]])]
    QAIC_values = current_QAIC_row[c(nested_QAIC_models_selected[[i]])]
    
    ordered_AIC_deltas = nested_AIC_models_selected[[i]][order(AIC_values, decreasing = FALSE)]
    ordered_QAIC_deltas = nested_QAIC_models_selected[[i]][order(QAIC_values, decreasing = FALSE)]
    
    
    final_AIC_list = ordered_AIC_deltas
    final_QAIC_list = ordered_QAIC_deltas
    
    for(j in 1:length(ordered_AIC_deltas))
    {
      AIC_removal_list = intersect(which(simplicity_key[ordered_AIC_deltas[j], ] == 1), ordered_AIC_deltas[j:length(ordered_AIC_deltas)])
      QAIC_removal_list = intersect(which(simplicity_key[ordered_QAIC_deltas[j], ] == 1), ordered_QAIC_deltas[j:length(ordered_QAIC_deltas)])
      
      
      # If a more complicated model has a greater AIC value:
      if(length(AIC_removal_list) > 0)
      {
        final_AIC_list = final_AIC_list[!final_AIC_list %in% AIC_removal_list]
      }
      
      if(length(QAIC_removal_list) > 0)
      {
        final_QAIC_list = final_QAIC_list[!final_QAIC_list %in% QAIC_removal_list]
      }
    }
    
    # Overwrite the original list with the updated list (nested models removed)
    nested_AIC_models_selected[[i]] = final_AIC_list
    nested_QAIC_models_selected[[i]] = final_QAIC_list
  }
  
  # Create an empty array to store the number of times each model is selected via the threshold selection rule. 
  nested_AIC_selection_times = rep(0, times = 8)
  delta_AIC_selection_times = rep(0, times = 8)
  nested_QAIC_selection_times = rep(0, times = 4)
  delta_QAIC_selection_times = rep(0, times = 4)
  
  for(i in 1:length(nested_AIC_models_selected))
  {
    # Increment the selection counts for the models selected in the i'th fit:
    nested_AIC_selection_times[nested_AIC_models_selected[[i]]] = nested_AIC_selection_times[nested_AIC_models_selected[[i]]] + 1
    delta_AIC_selection_times[delta_AIC_models_selected[[i]]] = delta_AIC_selection_times[delta_AIC_models_selected[[i]]] + 1
    
    nested_QAIC_selection_times[delta_QAIC_models_selected[[i]]] = nested_QAIC_selection_times[delta_QAIC_models_selected[[i]]] + 1
    delta_QAIC_selection_times[delta_QAIC_models_selected[[i]]] = delta_QAIC_selection_times[delta_QAIC_models_selected[[i]]] + 1
  }
  
  # Calculate the proportion of times each model was selected. 
  nested_AIC_selection_prop = nested_AIC_selection_times / num_reps
  delta_AIC_selection_prop = delta_AIC_selection_times / num_reps
  
  nested_QAIC_selection_prop = nested_QAIC_selection_times / num_reps
  delta_QAIC_selection_prop = delta_QAIC_selection_times / num_reps
  
  # AIC_table = data.frame(
  #   model = paste("M", seq(1:6), " (AIC)", sep = ""),
  #   delta_AIC = delta_AIC_selection_prop
  #   #nested_AIC = nested_AIC_selection_prop
  # )
  # 
  # QAIC_table = data.frame(
  #   model = paste("M", seq(1:3), " (QAIC)", sep = ""),
  #   delta_QAIC = delta_QAIC_selection_prop
  #   #nested_QAIC = nested_QAIC_selection_prop
  # )
  
  return_list = vector("list", length = 2) # Initialise a list to return from the function.
  
  return_list[[1]] = delta_AIC_selection_prop
  return_list[[2]] = delta_QAIC_selection_prop
  
  return(return_list)
}

threshold_sequence = seq(from = 0, to = 8, by = 0.5)

threshold_sequence = seq(from = 0, to = 8, by = 0.5)

AIC_table = data.frame(
  model = paste("(AIC) M", seq(1:8), sep = "")
)

QAIC_table = data.frame(
  model = paste("(QAIC) M", seq(1:4), sep = "")
)

for(i in 1:length(threshold_sequence))
{
  
  selection_output = prob_model_selected(AIC_threshold = threshold_sequence[i],
                                         QAIC_threshold = threshold_sequence[i],
                                         AIC_mat_delta = AIC_mat_delta,
                                         QAIC_mat_delta = QAIC_mat_delta)
  
  AIC_table[i + 1] = selection_output[[1]]
  QAIC_table[i + 1] = selection_output[[2]]
}

colnames(AIC_table)[2:ncol(AIC_table)] = threshold_sequence
colnames(QAIC_table)[2:ncol(QAIC_table)] = threshold_sequence

AIC_table_wide = gather(data = AIC_table, key = "Threshold", value = "Probability model is selected", 2:18)
QAIC_table_wide = gather(data = QAIC_table, key = "Threshold", value = "Probability model is selected", 2:18)

AIC_table_wide = AIC_table_wide[AIC_table_wide$model %in% c("(AIC) M3", "(AIC) M4", "(AIC) M7","(AIC) M8"), ]
QAIC_table_wide = QAIC_table_wide[QAIC_table_wide$model %in% c("(QAIC) M3", "(QAIC) M4"), ]


plot = ggplot() + 
  geom_point(data = AIC_table_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) + 
  geom_point(data = QAIC_table_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(0,15,1,19,1,19)) +
  guides(shape = guide_legend(override.aes = list(size=c(5,5,5,5,2,2))))

plot$layers[[1]]$aes_params$size = 5 # change the size of the geom_line layer from 1 to 0.5
plot$layers[[2]]$aes_params$size = 3   # change the size of the geom_point layer from 3 to 1
plot

full_plot = ggplot() + 
  geom_point(data = AIC_table_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) + 
  geom_point(data = QAIC_table_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(0,15,1,19,1,19)) +
  theme(legend.position="none")

full_plot$layers[[1]]$aes_params$size = 5 # change the size of the geom_line layer from 1 to 0.5
full_plot$layers[[2]]$aes_params$size = 3   # change the size of the geom_point layer from 3 to 1

plot_AIC = ggplot() + 
  geom_point(data = AIC_table_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(0,15,1,19), name = "AIC") +
  guides(shape = guide_legend(override.aes = list(size=c(5,5,5,5,5,5))))

plot_AIC$layers[[1]]$aes_params$size = 5 # change the size of the geom_line layer from 1 to 0.5


plot_QAIC = ggplot() + 
  geom_point(data = QAIC_table_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(1,19), name = "QAIC") +
  guides(shape = guide_legend(override.aes = list(size=c(2,2,2))))

plot_QAIC$layers[[1]]$aes_params$size = 3 # change the size of the geom_line layer from 1 to 0.5

plot_grid(
  full_plot
  , plot_grid(
    get_legend(plot_AIC)
    , get_legend(plot_QAIC)
    , nrow = 1
  )
  , nrow = 2
  , rel_heights = c(8,2)
)



# Calculate mean KLD and the associated standard errors. 
EKLD = data.frame(
  EKLD = colMeans(KLD_mat),
  Model = c(1:8),
  SE_EKLD = matrixStats::colSds(KLD_mat) / sqrt(nrow(KLD_mat))
)

EKLD[order(EKLD$EKLD), ]