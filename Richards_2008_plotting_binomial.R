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


# Constraint sets for later optimisation:
A_LL_M1 = matrix(c(1), 1, 1)
B_LL_M1 = matrix(c(0), 1, 1)

A_LL_M2 = matrix(c(1,0,0,1), 2,2)
B_LL_M2 = matrix(c(0,0), 2, 1)

A_LL_M3 = matrix(c(1,0,0,1), 2,2)
B_LL_M3 = matrix(c(0,0), 2, 1)

A_LL_M4 = matrix(c(1,0,0,1), 2,2)
B_LL_M4 = matrix(c(0,0), 2, 1)

A_LL_M5 = matrix(c(1,0,0,0,1,0,0,0,1), 3,3)
B_LL_M5 = matrix(c(0,0, 0), 3, 1)

A_LL_M6 = matrix(c(1,0,0,0,1,0,0,0,1), 3,3)
B_LL_M6 = matrix(c(0 , 0, 0), 3, 1) # Restrict 6 on 0. 

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
    
    number_AIC_models_selected[i] = length(nested_AIC_models_selected[[i]])
    number_nested_models_selected[i] = length(nested_AIC_models_selected[[i]])
    
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
  
  return_list = vector("list", length = 4) # Initialise a list to return from the function.
  
  return_list[[1]] = delta_AIC_selection_prop
  return_list[[2]] = delta_QAIC_selection_prop
  
  return_list[[3]] = number_AIC_models_selected
  return_list[[4]] = number_nested_models_selected
  
  return(return_list)
}


AIC_A = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_a.xlsx", sheetName = "AIC", header = FALSE))
KLD_estimates_A = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_a.xlsx", sheetName = "KLD_estimates", header = FALSE))
AIC_KLD_estimates_A = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_a.xlsx", sheetName = "AIC_KLD_estimates", header = FALSE))
QAIC_A = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_a.xlsx", sheetName = "QAIC", header = FALSE))

AIC_B = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_b.xlsx", sheetName = "AIC", header = FALSE))
KLD_estimates_B = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_b.xlsx", sheetName = "KLD_estimates", header = FALSE))
AIC_KLD_estimates_B = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_b.xlsx", sheetName = "AIC_KLD_estimates", header = FALSE))
QAIC_B = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_b.xlsx", sheetName = "QAIC", header = FALSE))

AIC_C = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_c.xlsx", sheetName = "AIC", header = FALSE))
KLD_estimates_C = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_c.xlsx", sheetName = "KLD_estimates", header = FALSE))
AIC_KLD_estimates_C = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_c.xlsx", sheetName = "AIC_KLD_estimates", header = FALSE))
QAIC_C = as.matrix(xlsx::read.xlsx(file = "binomial_scenario_c.xlsx", sheetName = "QAIC", header = FALSE))

AIC_A_delta = AIC_A - rowMins(AIC_A)
AIC_B_delta = AIC_B - rowMins(AIC_B)
AIC_C_delta = AIC_C - rowMins(AIC_C)


QAIC_A = QAIC_A[, 1:3]
QAIC_B = QAIC_B[, 1:3]
QAIC_C = QAIC_C[, 1:3]

QAIC_A_delta = QAIC_A - rowMins(QAIC_A)
QAIC_B_delta = QAIC_B - rowMins(QAIC_B)
QAIC_C_delta = QAIC_C - rowMins(QAIC_C)


threshold_sequence = seq(from = 0, to = 8, by = 0.5)

AIC_table_A = data.frame(
  model = paste("(AIC) M", seq(1:6), sep = "")
)
QAIC_table_A = data.frame(
  model = paste("(QAIC) M", seq(1:3), sep = "")
)

AIC_table_B = data.frame(
  model = paste("(AIC) M", seq(1:6), sep = "")
)
QAIC_table_B = data.frame(
  model = paste("(QAIC) M", seq(1:3), sep = "")
)

AIC_table_C = data.frame(
  model = paste("(AIC) M", seq(1:6), sep = "")
)
QAIC_table_C = data.frame(
  model = paste("(QAIC) M", seq(1:3), sep = "")
)

for(i in 1:length(threshold_sequence))
{
  # Scenario A:
  selection_output_A = prob_model_selected(AIC_threshold = threshold_sequence[i],
                                         QAIC_threshold = threshold_sequence[i],
                                         AIC_mat_delta = AIC_A_delta,
                                         QAIC_mat_delta = QAIC_A_delta)
  
  AIC_table_A[i + 1] = selection_output_A[[1]]
  QAIC_table_A[i + 1] = selection_output_A[[2]]
  
  # Scenario B: 
  selection_output_B = prob_model_selected(AIC_threshold = threshold_sequence[i],
                                         QAIC_threshold = threshold_sequence[i],
                                         AIC_mat_delta = AIC_B_delta,
                                         QAIC_mat_delta = QAIC_B_delta)
  
  AIC_table_B[i + 1] = selection_output_B[[1]]
  QAIC_table_B[i + 1] = selection_output_B[[2]]
  
  # Scenario C:
  selection_output_C = prob_model_selected(AIC_threshold = threshold_sequence[i],
                                         QAIC_threshold = threshold_sequence[i],
                                         AIC_mat_delta = AIC_C_delta,
                                         QAIC_mat_delta = QAIC_C_delta)
  
  AIC_table_C[i + 1] = selection_output_C[[1]]
  QAIC_table_C[i + 1] = selection_output_C[[2]]
}


selection_output_A = prob_model_selected(AIC_threshold = 6,
                                         QAIC_threshold = 6,
                                         AIC_mat_delta = AIC_A_delta,
                                         QAIC_mat_delta = QAIC_A_delta)


# Output the average number of models selected 
mean(selection_output_A[[3]])
mean(selection_output_A[[4]])




colnames(AIC_table_A)[2:ncol(AIC_table_A)] = threshold_sequence
colnames(QAIC_table_A)[2:ncol(QAIC_table_A)] = threshold_sequence

colnames(AIC_table_B)[2:ncol(AIC_table_B)] = threshold_sequence
colnames(QAIC_table_B)[2:ncol(QAIC_table_B)] = threshold_sequence

colnames(AIC_table_C)[2:ncol(AIC_table_C)] = threshold_sequence
colnames(QAIC_table_C)[2:ncol(QAIC_table_C)] = threshold_sequence



AIC_table_A_wide = gather(data = AIC_table_A, key = "Threshold", value = "Probability model is selected", 2:18)
QAIC_table_A_wide = gather(data = QAIC_table_A, key = "Threshold", value = "Probability model is selected", 2:18)

AIC_table_B_wide = gather(data = AIC_table_B, key = "Threshold", value = "Probability model is selected", 2:18)
QAIC_table_B_wide = gather(data = QAIC_table_B, key = "Threshold", value = "Probability model is selected", 2:18)

AIC_table_C_wide = gather(data = AIC_table_C, key = "Threshold", value = "Probability model is selected", 2:18)
QAIC_table_C_wide = gather(data = QAIC_table_C, key = "Threshold", value = "Probability model is selected", 2:18)

plot = ggplot() + 
  geom_point(data = AIC_table_A_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) + 
  geom_point(data = QAIC_table_A_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(0,7,15,1,13,19,1,13,19)) +
  guides(shape = guide_legend(override.aes = list(size=c(5,5,5,5,5,5,2,2,2))))

plot$layers[[1]]$aes_params$size = 5 # change the size of the geom_line layer from 1 to 0.5
plot$layers[[2]]$aes_params$size = 3   # change the size of the geom_point layer from 3 to 1
plot



plot_AIC_A = ggplot() + 
  geom_point(data = AIC_table_A_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(0,7,15,1,13,19), name = "AIC") +
  guides(shape = guide_legend(override.aes = list(size=c(5,5,5,5,5,5)), nrow = 3))

plot_AIC_A$layers[[1]]$aes_params$size = 5 # change the size of the geom_line layer from 1 to 0.5


plot_QAIC_A = ggplot() + 
  geom_point(data = QAIC_table_A_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(1,13,19), name = "QAIC") +
  guides(shape = guide_legend(override.aes = list(size=c(2,2,2))))

plot_QAIC_A$layers[[1]]$aes_params$size = 3 # change the size of the geom_line layer from 1 to 0.5



full_plot_A = ggplot() + 
  geom_point(data = AIC_table_A_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) + 
  geom_point(data = QAIC_table_A_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(0,7,15,1,13,19,1,13,19)) +
  xlab("AIC or QAIC Threshold") +
  scale_x_discrete(breaks=c(0:8)) +
  theme_bw() +
  theme(legend.position="none")

full_plot_A$layers[[1]]$aes_params$size = 4 # change the size of the geom_line layer from 1 to 0.5
full_plot_A$layers[[2]]$aes_params$size = 2   # change the size of the geom_point layer from 3 to 1


full_plot_B = ggplot() + 
  geom_point(data = AIC_table_B_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) + 
  geom_point(data = QAIC_table_B_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(0,7,15,1,13,19,1,13,19)) +
  xlab("AIC or QAIC Threshold") +
  scale_x_discrete(breaks=c(0:8)) +
  ylab("") + 
  theme_bw() +
  theme(legend.position="none")

full_plot_B$layers[[1]]$aes_params$size = 4 # change the size of the geom_line layer from 1 to 0.5
full_plot_B$layers[[2]]$aes_params$size = 2   # change the size of the geom_point layer from 3 to 1


full_plot_C = ggplot() + 
  geom_point(data = AIC_table_C_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) + 
  geom_point(data = QAIC_table_C_wide, mapping = aes(x = Threshold, y = `Probability model is selected`, shape = model, size = 1)) +
  scale_shape_manual(values = c(0,7,15,1,13,19,1,13,19)) +
  xlab("AIC or QAIC Threshold") +
  scale_x_discrete(breaks=c(0:8)) +
  ylab("") + 
  theme_bw() +
  theme(legend.position="none")

full_plot_C$layers[[1]]$aes_params$size = 4 # change the size of the geom_line layer from 1 to 0.5
full_plot_C$layers[[2]]$aes_params$size = 2   # change the size of the geom_point layer from 3 to 1


full_plot_grid = plot_grid(
  full_plot_A,
  full_plot_B,
  full_plot_C,
  nrow = 1
)

bottom_plot = plot_grid(
  full_plot_grid,
  plot_grid(
    NULL,
    get_legend(plot_AIC_A),
    get_legend(plot_QAIC_A),
    NULL,
    nrow = 1
  )
  , nrow = 2
  , rel_heights = c(3,1)
)


# Calculate mean KLD and the associated standard errors. 
EKLD_A = data.frame(
  Model = paste("M", c(1:6), sep = ""),
  AIC = colMeans(AIC_A),
  SE_AIC = matrixStats::colSds(AIC_A) / sqrt(nrow(AIC_A)),
  SD_AIC = matrixStats:::colSds(AIC_A),
  
  EKLD_AIC = colMeans(AIC_KLD_estimates_A),
  SD_EKLD_AIC = matrixStats::colSds(AIC_KLD_estimates_A),
  SE_EKLD_AIC = matrixStats::colSds(AIC_KLD_estimates_A) / sqrt(nrow(AIC_KLD_estimates_A))
)

# Calculate mean KLD and the associated standard errors. 
EKLD_B = data.frame(
  Model = paste("M", c(1:6), sep = ""),
  AIC = colMeans(AIC_B),
  SE_AIC = matrixStats::colSds(AIC_B) / sqrt(nrow(AIC_B)),
  SD_AIC = matrixStats:::colSds(AIC_B),
  
  EKLD_AIC = colMeans(AIC_KLD_estimates_B),
  SD_EKLD_AIC = matrixStats::colSds(AIC_KLD_estimates_B),
  SE_EKLD_AIC = matrixStats::colSds(AIC_KLD_estimates_B) / sqrt(nrow(AIC_KLD_estimates_B))
)

# Calculate mean KLD and the associated standard errors. 
EKLD_C = data.frame(
  Model = paste("M", c(1:6), sep = ""),
  AIC = colMeans(AIC_C),
  SE_AIC = matrixStats::colSds(AIC_C) / sqrt(nrow(AIC_C)),
  SD_AIC = matrixStats:::colSds(AIC_C),
  
  EKLD_AIC = colMeans(AIC_KLD_estimates_C),
  SD_EKLD_AIC = matrixStats::colSds(AIC_KLD_estimates_C),
  SE_EKLD_AIC = matrixStats::colSds(AIC_KLD_estimates_C) / sqrt(nrow(AIC_KLD_estimates_C))
)

scen_A = ggplot(data = EKLD_A, mapping = aes(x = EKLD_AIC, y = AIC)) + 
  geom_point() + 
  geom_errorbar(
    aes(
      ymin = AIC - SD_AIC,
      ymax = AIC + SD_AIC),
    width = 0.2
  ) + 
  ylim(NA, max(EKLD_A$AIC + EKLD_A$SD_AIC + 5)) + 
  ggrepel::geom_text_repel(
    label = EKLD_A$Model, 
    nudge_x = 0,
    nudge_y = (EKLD_A$SD_AIC + 5),
    segment.color = NA
  ) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size = 0.4) + 
  xlab(latex2exp::TeX("$2(E_p\\[I(\\textbf{p},\\textbf{π})\\] - c)$")) + 
  ylab("AIC") + 
  ggtitle(label = "Scenario A \nStrong kelp effect \nHigh overdispersion") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) 

scen_B = ggplot(data = EKLD_B, mapping = aes(x = EKLD_AIC, y = AIC)) + 
  geom_point() + 
  geom_errorbar(
    aes(
      ymin = AIC - SD_AIC,
      ymax = AIC + SD_AIC),
    width = 0.2
  ) + 
  ylim(NA, max(EKLD_B$AIC + EKLD_B$SD_AIC + 4)) + 
  ggrepel::geom_text_repel(
    label = EKLD_B$Model, 
    nudge_x = 0,
    nudge_y = (EKLD_B$SD_AIC * 1.5),
    segment.color = NA
  ) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size = 0.4) + 
  xlab(latex2exp::TeX("$2(E_p\\[I(\\textbf{p},\\textbf{π})\\] - c)$")) + 
  ylab("AIC") + 
  ggtitle(label = "Scenario B \nWeak kelp effect \nHigh overdispersion") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) 


scen_C = ggplot(data = EKLD_C, mapping = aes(x = EKLD_AIC, y = AIC)) + 
  geom_point() + 
  geom_errorbar(
    aes(
      ymin = AIC - SD_AIC,
      ymax = AIC + SD_AIC),
    width = 0.2
  ) + 
  ylim(NA, max(EKLD_C$AIC + EKLD_C$SD_AIC + 5)) + 
  ggrepel::geom_text_repel(
    label = EKLD_C$Model, 
    nudge_x = 0,
    nudge_y = (EKLD_C$SD_AIC * 1.5),
    segment.color = NA
  ) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size = 0.4) + 
  xlab(latex2exp::TeX("$2(E_p\\[I(\\textbf{p},\\textbf{π})\\] - c)$")) + 
  ylab("AIC") + 
  ggtitle(label = "Scenario C \nStrong kelp effect \nLow overdispersion") + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5))


top_plot = plot_grid(
  scen_A,
  scen_B,
  scen_C,
  nrow = 1
)

final_plot = plot_grid(
  top_plot,
  bottom_plot,
  nrow = 2 
)
final_plot # Going to do a manual save on this one to avoid compression. 
