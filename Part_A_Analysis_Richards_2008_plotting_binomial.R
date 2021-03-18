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

# Binomial Example: 

kelp_bunches = seq(1:4) # The number of kelp bunches. 
tau = 15 # The duration of the experiment. 
alpha = 0.075 # Pre-defined capture rate (per hour). 
beta = 0.5 # The effectiveness of the kelp in mitigating perch capture. 
num_replicates = 5 # The number of truth values to be generated and which the models will be fit to (for each treatment). 
num_perch = 10 # The number of prey fish in the tank. 
phi = 0.025 # The pre-determined overdispersion factor.
num_reps = 1000 # The number of times the KLD will be calculated. 
Z = 10000 # The number of new points to use in calculation of cumulative IC (for each treatment). 
AIC_threshold = 6 # The threshold for delta_AIC-based selection. 
QAIC_threshold = 6


simplicity_key = matrix(
  data = c(
    c(0,1,1,1,1,1),
    c(0,0,0,0,1,1),
    c(0,0,0,0,1,1),
    c(0,0,0,0,1,1),
    c(0,0,0,0,0,0),
    c(0,0,0,0,0,0)
  ), nrow = 6, byrow = TRUE)



prob_model_selected = function(AIC_threshold, QAIC_threshold, AIC_mat_delta, QAIC_mat_delta)
{
  # Create a list of models selected from each fit.
  nested_AIC_models_selected = list() 
  delta_AIC_models_selected = list()
  delta_QAIC_models_selected = list()
  nested_QAIC_models_selected = list()
  
  number_AIC_models_selected = c()
  number_nested_models_selected = c()
  
  best_model_selected_delta = 0
  best_model_selected_nested = 0 
  
  for(i in 1:nrow(AIC_mat_delta))
  {
    # Extract the i'th rows. 
    current_AIC_row = AIC_mat_delta[i, ] 
    current_QAIC_row = QAIC_mat_delta[i, ]
    
    # Determine which models are below the delta-AIC threshold. 
    nested_AIC_models_selected[[i]] = which(current_AIC_row <= AIC_threshold) 
    delta_AIC_models_selected[[i]] = which(current_AIC_row <= AIC_threshold)
    delta_QAIC_models_selected[[i]] = which(current_QAIC_row <= QAIC_threshold)
    nested_QAIC_models_selected[[i]] = which(current_QAIC_row <= QAIC_threshold)
    
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
  nested_AIC_selection_times = rep(0, times = 6)
  delta_AIC_selection_times = rep(0, times = 6)
  nested_QAIC_selection_times = rep(0, times = 3)
  delta_QAIC_selection_times = rep(0, times = 3)
  
  for(i in 1:length(nested_AIC_models_selected))
  {
    # Increment the selection counts for the models selected in the i'th fit:
    nested_AIC_selection_times[nested_AIC_models_selected[[i]]] = nested_AIC_selection_times[nested_AIC_models_selected[[i]]] + 1
    delta_AIC_selection_times[delta_AIC_models_selected[[i]]] = delta_AIC_selection_times[delta_AIC_models_selected[[i]]] + 1
    
    nested_QAIC_selection_times[delta_QAIC_models_selected[[i]]] = nested_QAIC_selection_times[delta_QAIC_models_selected[[i]]] + 1
    delta_QAIC_selection_times[delta_QAIC_models_selected[[i]]] = delta_QAIC_selection_times[delta_QAIC_models_selected[[i]]] + 1
    
    number_AIC_models_selected[i] = length(delta_AIC_models_selected[[i]])
    number_nested_models_selected[i] = length(nested_AIC_models_selected[[i]])
    
    
    if(5 %in% delta_AIC_models_selected[[i]])
    {
      best_model_selected_delta = best_model_selected_delta + 1
    }
    
    if(5 %in% nested_AIC_models_selected[[i]])
    {
      best_model_selected_nested = best_model_selected_nested + 1
    }
    
  }
  
  # Calculate the proportion of times each model was selected. 
  nested_AIC_selection_prop = nested_AIC_selection_times / num_reps
  delta_AIC_selection_prop = delta_AIC_selection_times / num_reps
  
  nested_QAIC_selection_prop = nested_QAIC_selection_times / num_reps
  delta_QAIC_selection_prop = delta_QAIC_selection_times / num_reps
  
  return_list = vector("list", length = 6) # Initialise a list to return from the function.
  
  return_list[[1]] = delta_AIC_selection_prop
  return_list[[2]] = delta_QAIC_selection_prop
  
  return_list[[3]] = number_AIC_models_selected
  return_list[[4]] = number_nested_models_selected
  
  return_list[[5]] = best_model_selected_delta / nrow(AIC_mat_delta)
  return_list[[6]] = best_model_selected_nested / nrow(AIC_mat_delta)
  
  return(return_list)
}

# Read in the scenario data: 
AIC_A = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_a.xlsx", sheetName = "AIC", header = FALSE))
KLD_estimates_A = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_a.xlsx", sheetName = "KLD_estimates", header = FALSE))
AIC_KLD_estimates_A = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_a.xlsx", sheetName = "AIC_KLD_estimates", header = FALSE))
QAIC_A = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_a.xlsx", sheetName = "QAIC", header = FALSE))


AIC_A_delta = AIC_A - rowMins(AIC_A)
QAIC_A = QAIC_A[, 1:3]
QAIC_A_delta = QAIC_A - rowMins(QAIC_A)

threshold_sequence = seq(from = 0, to = 8, by = 0.5)

AIC_table_A = data.frame(
  model = paste("(AIC) M", seq(1:6), sep = "")
)
QAIC_table_A = data.frame(
  model = paste("(QAIC) M", seq(1:3), sep = "")
)


for(i in 1:length(threshold_sequence))
{
  # Scenario A:
  selection_output = prob_model_selected(AIC_threshold = threshold_sequence[i],
                                         QAIC_threshold = threshold_sequence[i],
                                         AIC_mat_delta = AIC_A_delta,
                                         QAIC_mat_delta = QAIC_A_delta)
  
  AIC_table_A[i + 1] = selection_output[[1]]
  QAIC_table_A[i + 1] = selection_output[[2]]
}


selection_output = prob_model_selected(AIC_threshold = 6,
                                         QAIC_threshold = 6,
                                         AIC_mat_delta = AIC_A_delta,
                                         QAIC_mat_delta = QAIC_A_delta)


# Output the average number of models selected 
mean(selection_output[[3]])
mean(selection_output[[4]])

# Output the EKLD estimate for models 5 and 6. 
colmeans(KLD_estimates_A[, 5:6])

selection_output[[5]]
selection_output[[6]]