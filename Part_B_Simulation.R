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

predator_1 = +0.1 # The reduction in effectiveness associated with the second predator. 
phenotype_1 = -0.05 # The reduction in effectiveness associated with the second phenotype.  


num_replicates = 5 # The number of truth values to be generated and which the models will be fit to (for each treatment). 
num_perch = 20 # The number of prey fish in the tank. 
phi = 0.025 # The pre-determined overdispersion factor.
num_reps = 1000 # The number of times the KLD will be calculated. 
Z = 10000 # The number of new points to use in calculation of cumulative IC (for each treatment). 
AIC_threshold = 6 # The threshold for delta_AIC-based selection. 
QAIC_threshold = 6


# Constraint sets for later optimisation:
A_LL_M1 = matrix(c(1), 1, 1)
B_LL_M1 = matrix(c(0), 1, 1)

A_LL_M2_1 = matrix(c(1,0,0,1), 2,2)
B_LL_M2_1 = matrix(c(0,0), 2, 1)

A_LL_M2_2_SA_PA = matrix(c(1,0,0,1), 2,2)
B_LL_M2_2_SA_PA = matrix(c(0,.Machine$integer.max), 2,1)

A_LL_M2_2_SA_PB = matrix(c(1,0,0,1), 2,2)
B_LL_M2_2_SA_PB = matrix(c(0,.Machine$integer.max), 2,1)

A_LL_M2_2_SB_PA = matrix(c(1,0,0,1), 2,2)
B_LL_M2_2_SB_PA = matrix(c(0,.Machine$integer.max), 2,1)

A_LL_M2_2_SB_PB = matrix(c(1,0,0,1), 2,2)
B_LL_M2_2_SB_PB = matrix(c(0,.Machine$integer.max), 2,1)

A_LL_M2_2_PRED = matrix(c(1,0,0,1), 2,2)
B_LL_M2_2_PRED = matrix(c(0,.Machine$integer.max), 2,1)

A_LL_M2_2_PHEN = matrix(c(1,0,0,1), 2,2)
B_LL_M2_2_PHEN = matrix(c(0,.Machine$integer.max), 2,1)

A_LL_M3_1 = matrix(c(1,0,0,1), 2,2)
B_LL_M3_1 = matrix(c(0,0), 2, 1)

A_LL_M3_2 = matrix(c(1,0,0,1), 2,2)
B_LL_M3_2 = matrix(c(0,.Machine$integer.max), 2, 1)

A_LL_M3_3 = matrix(c(1,0,0,0,1,0,0,0,1), 3,3)
B_LL_M3_3 = matrix(c(0, 0, .Machine$integer.max), 3, 1)

A_LL_M3_4 = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), 4,4)
B_LL_M3_4 = matrix(c(0, 0, .Machine$integer.max, .Machine$integer.max ), 4, 1)

A_LL_M3_5 = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1), 4,4)
B_LL_M3_5 = matrix(c(0, 0, .Machine$integer.max, .Machine$integer.max ), 4, 1)

A_LL_M3_6 = matrix(c(
  c(1,0,0,0,0,0,0,0),
  c(0,1,0,0,0,0,0,0),
  c(0,0,1,0,0,0,0,0),
  c(0,0,0,1,0,0,0,0),
  c(0,0,0,0,1,0,0,0),
  c(0,0,0,0,0,1,0,0),
  c(0,0,0,0,0,0,1,0),
  c(0,0,0,0,0,0,0,1)
), 8,8)
B_LL_M3_6 = matrix(c(0, 0, .Machine$integer.max, .Machine$integer.max, .Machine$integer.max,
                     .Machine$integer.max, .Machine$integer.max, .Machine$integer.max), 8, 1)

# CONTRAINTS FOR BETA-BINOMIAL MODELS:

A_LL_M4 = diag(2)
B_LL_M4 = matrix(c(0,0), 2, 1)

A_LL_M5_1 = diag(3)
B_LL_M5_1 = matrix(c(0,0,0), 3, 1)

A_LL_M5_2_SA_PA = diag(3)
B_LL_M5_2_SA_PA = matrix(c(0,.Machine$integer.max,0), 3,1)

A_LL_M5_2_SA_PB = diag(3)
B_LL_M5_2_SA_PB = matrix(c(0,.Machine$integer.max,0), 3,1)

A_LL_M5_2_SB_PA = diag(3)
B_LL_M5_2_SB_PA = matrix(c(0,.Machine$integer.max,0), 3,1)

A_LL_M5_2_SB_PB = diag(3)
B_LL_M5_2_SB_PB = matrix(c(0,.Machine$integer.max,0), 3,1)

A_LL_M5_2_PRED = diag(3)
B_LL_M5_2_PRED = matrix(c(0,.Machine$integer.max,0), 3,1)

A_LL_M5_2_PHEN = diag(3)
B_LL_M5_2_PHEN = matrix(c(0,.Machine$integer.max,0), 3,1)

A_LL_M6_1 = diag(3)
B_LL_M6_1 = matrix(c(0,0,0), 3, 1)

A_LL_M6_2 = diag(3)
B_LL_M6_2 = matrix(c(0,.Machine$integer.max, 0), 3, 1)

A_LL_M6_3 = diag(4)
B_LL_M6_3 = matrix(c(0, 0, .Machine$integer.max, 0), 4, 1)

A_LL_M6_4 = diag(5)
B_LL_M6_4 = matrix(c(0, 0, .Machine$integer.max, .Machine$integer.max, 0 ), 5, 1)

A_LL_M6_5 = diag(5)
B_LL_M6_5 = matrix(c(0, 0, .Machine$integer.max, .Machine$integer.max, 0), 5, 1)

A_LL_M6_6 = diag(9)
B_LL_M6_6 = matrix(c(0, 0, .Machine$integer.max, .Machine$integer.max, .Machine$integer.max,
                     .Machine$integer.max, .Machine$integer.max, .Machine$integer.max, 0), 9, 1)







simplicity_key = matrix(
  data = c(
    c(0,1,1,1,1,1),
    c(0,0,0,0,1,1),
    c(0,0,0,0,1,1),
    c(0,0,0,0,1,1),
    c(0,0,0,0,0,0),
    c(0,0,0,0,0,0)
  ), nrow = 6, byrow = TRUE)


# Collection of functions associated with each of the mean probabilities of survival 
# for each model.
# @param x: The number of bunches in the pool (in each relevant function, below).

# Generate a list of x-values corresponding to the different treatments: 
x = rep(1:4, each = num_replicates, times = 4)
phenotype = rep(0:1, each = num_replicates * 8)
predator = rep(0:1, each = num_replicates * 4, times = 2)

# Bind each of the variables into one matrix:
variable_matrix = as.matrix(cbind(x,phenotype,predator))


# The true function - specified as LL_M3_5 in the log likelihood functions. 
p3 = function(kelp, phenotype, predator) {
  
  numerator = exp(-tau * alpha + beta * kelp + predator_1 * predator + phenotype_1 * predator)
  denominator = 1 + exp(-tau * alpha + beta * kelp + predator_1 * predator + phenotype_1 * predator)
  
  p_bar = numerator / denominator
  
  return(p_bar)
}

# We remove the log of the combination choose(num_perch, y) within each of M1 - M3
# as it is constant with respect to alpha in the maximisation of the log likelihood. 
LL_M1 = function(par, y_perch) {
  
  p_bar = exp(-tau * par[1])
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}

# Base case using only number of kelp:
LL_M2_1 = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}


# Base case using only number of kelp:
LL_M2_2_SA_PA = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}

# Base case using only number of kelp:
LL_M2_2_SB_PA = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}

# Base case using only number of kelp:
LL_M2_2_SA_PB = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}

# Base case using only number of kelp:
LL_M2_2_SB_PB = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}


# Base case using only predator type:
LL_M2_2_PRED = function(par, y_perch, x_predator) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_predator))
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}

# Base case using only phenotype:
LL_M2_2_PHEN = function(par, y_perch, x_phenotype) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_phenotype))
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}



# Base logistic case using only number of kelp:
LL_M3_1 = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  
  numerator = exp(-tau * alpha_LL + beta_LL * x_kelp)
  denominator = 1 + exp(-alpha_LL * tau + beta_LL * x_kelp)
  
  p_bar = numerator / denominator
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}


# Alternative base case using only the predator type: 
LL_M3_2 = function(par, y_perch, x_predator) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_predator)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_predator)
  
  p_bar = numerator / denominator
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}


# Case using number of kelp and phenotype:
LL_M3_3 = function(par, y_perch, x_kelp, x_phenotype) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  beta_LL_2 = par[3]
  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype)
  
  p_bar = numerator / denominator
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}


# Case using interaction between number of kelp and predator type:
LL_M3_4 = function(par, y_perch, x_kelp, x_predator) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  beta_LL_2 = par[3]
  beta_LL_3 = par[4]
  
  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_predator + beta_LL_3 * x_predator * x_kelp)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_predator + beta_LL_3 * x_predator * x_kelp)
  
  p_bar = numerator / denominator
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}


# Case using number of kelp, phenotype and predator type: 
LL_M3_5 = function(par, y_perch, x_kelp, x_phenotype, x_predator) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  beta_LL_2 = par[3]
  beta_LL_3 = par[4]

  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype + beta_LL_3 * x_predator)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype + beta_LL_3 * x_predator)
  
  p_bar = numerator / denominator
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}

# Case using number of kelp, phenotype and predator type, and full pair-wise 
# interactions and three-way interaction. This model represents an over-specification 
# of the model.
LL_M3_6 = function(par, y_perch, x_kelp, x_phenotype, x_predator) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  beta_LL_2 = par[3]
  beta_LL_3 = par[4]
  beta_LL_4 = par[5]
  beta_LL_5 = par[6]
  beta_LL_6 = par[7]
  beta_LL_7 = par[8]
  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype + 
                    beta_LL_3 * x_predator + beta_LL_4 * x_kelp * x_phenotype + 
                    beta_LL_5 * x_kelp * x_phenotype + beta_LL_6 * x_phenotype * x_kelp + 
                    beta_LL_7 * x_kelp * x_phenotype * x_predator)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype + 
                          beta_LL_3 * x_predator + beta_LL_4 * x_kelp * x_phenotype + 
                          beta_LL_5 * x_kelp * x_phenotype + beta_LL_6 * x_phenotype * x_kelp + 
                          beta_LL_7 * x_kelp * x_phenotype * x_predator)
  
  p_bar = numerator / denominator
  
  return(sum(log(choose(num_perch, y_perch)) + y_perch * log(p_bar) + (num_perch - y_perch) * log(1 - p_bar)))
}




# BETA-BINOMIAL FUNCTIONS:

LL_M4 = function(par, y_perch) {
  
  phi = par[2]
  p_bar = exp(-tau * par[1])
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


# Base case using only number of kelp:
LL_M5_1 = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  phi = par[3]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


# Base case using only number of kelp:
LL_M5_2_SA_PA = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  phi = par[3]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}

# Base case using only number of kelp:
LL_M5_2_SB_PA = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  phi = par[3]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}

# Base case using only number of kelp:
LL_M5_2_SA_PB = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  phi = par[3]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}

# Base case using only number of kelp:
LL_M5_2_SB_PB = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  phi = par[3]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_kelp))

  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


# Base case using only predator type:
LL_M5_2_PRED = function(par, y_perch, x_predator) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  phi = par[3]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_predator))
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}

# Base case using only phenotype:
LL_M5_2_PHEN = function(par, y_perch, x_phenotype) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  phi = par[3]
  
  p_bar = exp(-tau * alpha_LL / (1 + beta_LL * x_phenotype))
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}



# Base logistic case using only number of kelp:
LL_M6_1 = function(par, y_perch, x_kelp) {
  
  alpha_LL = par[1]
  beta_LL = par[2]
  phi = par[3]
  
  numerator = exp(-tau * alpha_LL + beta_LL * x_kelp)
  denominator = 1 + exp(-alpha_LL * tau + beta_LL * x_kelp)
  
  p_bar = numerator / denominator
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


# Alternative base case using only the predator type: 
LL_M6_2 = function(par, y_perch, x_predator) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  phi = par[3]
  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_predator)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_predator)
  
  p_bar = numerator / denominator
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


# Case using number of kelp and phenotype:
LL_M6_3 = function(par, y_perch, x_kelp, x_phenotype) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  beta_LL_2 = par[3]
  phi = par[4]
  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype)
  
  p_bar = numerator / denominator
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


# Case using interaction between number of kelp and predator type:
LL_M6_4 = function(par, y_perch, x_kelp, x_predator) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  beta_LL_2 = par[3]
  beta_LL_3 = par[4]
  phi = par[5]
  
  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_predator + beta_LL_3 * x_predator * x_kelp)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_predator + beta_LL_3 * x_predator * x_kelp)
  
  p_bar = numerator / denominator
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


# Case using number of kelp, phenotype and predator type: 
LL_M6_5 = function(par, y_perch, x_kelp, x_phenotype, x_predator) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  beta_LL_2 = par[3]
  beta_LL_3 = par[4]
  phi = par[5]
  
  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype + beta_LL_3 * x_predator)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype + beta_LL_3 * x_predator)
  
  p_bar = numerator / denominator
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}

# Case using number of kelp, phenotype and predator type, and full pair-wise 
# interactions and three-way interaction. This model represents an over-specification 
# of the model.
LL_M6_6 = function(par, y_perch, x_kelp, x_phenotype, x_predator) {
  
  alpha_LL = par[1]
  beta_LL_1 = par[2]
  beta_LL_2 = par[3]
  beta_LL_3 = par[4]
  beta_LL_4 = par[5]
  beta_LL_5 = par[6]
  beta_LL_6 = par[7]
  beta_LL_7 = par[8]
  phi = par[9]
  
  numerator = exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype + 
                    beta_LL_3 * x_predator + beta_LL_4 * x_kelp * x_phenotype + 
                    beta_LL_5 * x_kelp * x_phenotype + beta_LL_6 * x_phenotype * x_kelp + 
                    beta_LL_7 * x_kelp * x_phenotype * x_predator)
  denominator = 1 + exp(-tau * alpha_LL + beta_LL_1 * x_kelp + beta_LL_2 * x_phenotype + 
                          beta_LL_3 * x_predator + beta_LL_4 * x_kelp * x_phenotype + 
                          beta_LL_5 * x_kelp * x_phenotype + beta_LL_6 * x_phenotype * x_kelp + 
                          beta_LL_7 * x_kelp * x_phenotype * x_predator)
  
  p_bar = numerator / denominator
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y_perch + a) + lgamma(num_perch - y_perch + b) - 
               lgamma(y_perch +  1) - lgamma(num_perch - y_perch + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}





process_rep = function(i)
{
  
  # Generate the i'th "truth" dataset: 
  truth_set = c()
  
  for(j in 1:length(x))
  {
    p_bar = p3(kelp = x[j], phenotype = phenotype[j], predator = predator[j])
    
    alpha_truth = p_bar / phi
    beta_truth = (1 - p_bar) / phi
    
    truth_set[j] = extraDistr::rbbinom(n = 1, size = num_perch, 
                                       alpha = alpha_truth, beta = beta_truth)
    
  }
  
  # Fit each of the models to the generated data:
  
  M1_fit = maxLik::maxLik(logLik = LL_M1, start = c(1), constraints = list(ineqA = A_LL_M1, ineqB = B_LL_M1), method = "BFGS", y_perch = truth_set)
  
  M2_1_fit = maxLik::maxLik(logLik = LL_M2_1, start = c(1,1), constraints = list(ineqA = A_LL_M2_1, ineqB = B_LL_M2_1), y_perch = truth_set, x_kelp = x)
  
  M2_2_SA_PA_fit = maxLik::maxLik(logLik = LL_M2_2_SA_PA, start = c(1,1), constraints = list(ineqA = A_LL_M2_2_SA_PA, ineqB = B_LL_M2_2_SA_PA), 
                                  y_perch = truth_set[1:20], x_kelp = x[1:20])
    
  M2_2_SA_PB_fit = maxLik::maxLik(logLik = LL_M2_2_SA_PB, start = c(1,1), constraints = list(ineqA = A_LL_M2_2_SA_PB, ineqB = B_LL_M2_2_SA_PB), 
                                  y_perch = truth_set[21:40], x_kelp = x[21:40])
  
  M2_2_SB_PA_fit = maxLik::maxLik(logLik = LL_M2_2_SB_PA, start = c(1,1), constraints = list(ineqA = A_LL_M2_2_SB_PA, ineqB = B_LL_M2_2_SB_PA), 
                                  y_perch = truth_set[41:60], x_kelp = x[41:60])
  
  M2_2_SB_PB_fit = maxLik::maxLik(logLik = LL_M2_2_SB_PB, start = c(1,1), constraints = list(ineqA = A_LL_M2_2_SB_PB, ineqB = B_LL_M2_2_SB_PB), 
                                  y_perch = truth_set[61:80], x_kelp = x[61:80])
  
  M2_2_SA_PA_fit$maximum + M2_2_SA_PB_fit$maximum + M2_2_SB_PA_fit$maximum + M2_2_SB_PB_fit$maximum
  
  M2_2_PRED_fit = maxLik::maxLik(logLik = LL_M2_2_PRED, start = c(1,1), constraints = list(ineqA = A_LL_M2_2_PRED, ineqB = B_LL_M2_2_PRED), 
                                 y_perch = truth_set, x_predator = predator)
  
  M2_2_PHEN_fit = maxLik::maxLik(logLik = LL_M2_2_PHEN, start = c(1,1), constraints = list(ineqA = A_LL_M2_2_PHEN, ineqB = B_LL_M2_2_PHEN), 
                                 y_perch = truth_set, x_phenotype = phenotype)
  
  
  M3_1_fit = maxLik::maxLik(logLik = LL_M3_1, start = c(1,1), constraints = list(ineqA = A_LL_M3_1, ineqB = B_LL_M3_1), x_kelp = x, y_perch = truth_set)
  M3_2_fit = maxLik::maxLik(logLik = LL_M3_2, start = c(1,1), constraints = list(ineqA = A_LL_M3_2, ineqB = B_LL_M3_2), x_predator = predator, y_perch = truth_set)
  M3_3_fit = maxLik::maxLik(logLik = LL_M3_3, start = c(1,1,1), constraints = list(ineqA = A_LL_M3_3, ineqB = B_LL_M3_3), x_kelp = x, x_phenotype = phenotype, y_perch = truth_set)
  M3_4_fit = maxLik::maxLik(logLik = LL_M3_4, start = c(1,1,1,1), constraints = list(ineqA = A_LL_M3_4, ineqB = B_LL_M3_4), x_kelp = x, x_predator = predator, y_perch = truth_set, iterlim = 1000)
  M3_5_fit = maxLik::maxLik(logLik = LL_M3_5, start = c(1,1,1,1), constraints = list(ineqA = A_LL_M3_5, ineqB = B_LL_M3_5), y_perch = truth_set, x_kelp = x, x_phenotype = phenotype,  x_predator = predator, iterlim = 1000)
  M3_6_fit = maxLik::maxLik(logLik = LL_M3_6, start = c(1,1,1,1,1,1,1,1), constraints = list(ineqA = A_LL_M3_6, ineqB = B_LL_M3_6), y_perch = truth_set, x_kelp = x, x_phenotype = phenotype,  x_predator = predator, iterlim = 1000)
  
  # FIT THE BETA-BINOMIAL MODELS:
  M4_fit = maxLik::maxLik(logLik = LL_M4, start = c(1,1), constraints = list(ineqA = A_LL_M4, ineqB = B_LL_M4), method = "BFGS", y_perch = truth_set)
  
  M5_1_fit = maxLik::maxLik(logLik = LL_M5_1, start = c(1,1,1), constraints = list(ineqA = A_LL_M5_1, ineqB = B_LL_M5_1), y_perch = truth_set, x_kelp = x)
  
  M5_2_SA_PA_fit = maxLik::maxLik(logLik = LL_M5_2_SA_PA, start = c(1,1,1), constraints = list(ineqA = A_LL_M5_2_SA_PA, ineqB = B_LL_M5_2_SA_PA), 
                                  y_perch = truth_set[1:20], x_kelp = x[1:20])
  
  M5_2_SA_PB_fit = maxLik::maxLik(logLik = LL_M5_2_SA_PB, start = c(1,1,1), constraints = list(ineqA = A_LL_M5_2_SA_PB, ineqB = B_LL_M5_2_SA_PB), 
                                  y_perch = truth_set[21:40], x_kelp = x[21:40])
  
  M5_2_SB_PA_fit = maxLik::maxLik(logLik = LL_M5_2_SB_PA, start = c(1,1,1), constraints = list(ineqA = A_LL_M5_2_SB_PA, ineqB = B_LL_M5_2_SB_PA), 
                                  y_perch = truth_set[41:60], x_kelp = x[41:60])
  
  M5_2_SB_PB_fit = maxLik::maxLik(logLik = LL_M5_2_SB_PB, start = c(1,1,1), constraints = list(ineqA = A_LL_M5_2_SB_PB, ineqB = B_LL_M5_2_SB_PB), 
                                  y_perch = truth_set[61:80], x_kelp = x[61:80])
  
  # M5_2_SA_PA_fit$maximum + M5_2_SA_PB_fit$maximum + M5_2_SB_PA_fit$maximum + M5_2_SB_PB_fit$maximum
  
  M5_2_PRED_fit = maxLik::maxLik(logLik = LL_M5_2_PRED, start = c(1,1,1), constraints = list(ineqA = A_LL_M5_2_PRED, ineqB = B_LL_M5_2_PRED), 
                                 y_perch = truth_set, x_predator = predator)
  
  M5_2_PHEN_fit = maxLik::maxLik(logLik = LL_M5_2_PHEN, start = c(1,1,1), constraints = list(ineqA = A_LL_M5_2_PHEN, ineqB = B_LL_M5_2_PHEN), 
                                 y_perch = truth_set, x_phenotype = phenotype)
  
  
  M6_1_fit = maxLik::maxLik(logLik = LL_M6_1, start = c(1,1,1), constraints = list(ineqA = A_LL_M6_1, ineqB = B_LL_M6_1), x_kelp = x, y_perch = truth_set)
  M6_2_fit = maxLik::maxLik(logLik = LL_M6_2, start = c(1,1,1), constraints = list(ineqA = A_LL_M6_2, ineqB = B_LL_M6_2), x_predator = predator, y_perch = truth_set)
  M6_3_fit = maxLik::maxLik(logLik = LL_M6_3, start = c(1,1,1,1), constraints = list(ineqA = A_LL_M6_3, ineqB = B_LL_M6_3), x_kelp = x, x_phenotype = phenotype, y_perch = truth_set)
  M6_4_fit = maxLik::maxLik(logLik = LL_M6_4, start = c(1,1,1,1,1), constraints = list(ineqA = A_LL_M6_4, ineqB = B_LL_M6_4), x_kelp = x, x_predator = predator, y_perch = truth_set, iterlim = 1000)
  M6_5_fit = maxLik::maxLik(logLik = LL_M6_5, start = c(1,1,1,1,1), constraints = list(ineqA = A_LL_M6_5, ineqB = B_LL_M6_5), y_perch = truth_set, x_kelp = x, x_phenotype = phenotype,  x_predator = predator, iterlim = 1000)
  M6_6_fit = maxLik::maxLik(logLik = LL_M6_6, start = c(1,1,1,1,1,1,1,1,1), constraints = list(ineqA = A_LL_M6_6, ineqB = B_LL_M6_6), y_perch = truth_set, x_kelp = x, x_phenotype = phenotype,  x_predator = predator, iterlim = 1000)
  
  
  
  
  
  AIC_results = rep(0, times = 22)
  
  AIC_results[1] = -2 * M1_fit$maximum + 2 * length(M1_fit$estimate)
  AIC_results[2] = -2 * M2_1_fit$maximum + 2 * length(M2_1_fit$estimate) 
  AIC_results[3] = -2 * (M2_2_SA_PA_fit$maximum + M2_2_SA_PB_fit$maximum + 
                           M2_2_SB_PA_fit$maximum + M2_2_SB_PB_fit$maximum) + 
                           2 * 4 * length(M2_2_SA_PA_fit$estimate)
  
  AIC_results[4] = -2 * M2_2_PRED_fit$maximum + 2 * length(M2_2_PRED_fit$estimate)
  AIC_results[5] = -2 * M2_2_PHEN_fit$maximum + 2 * length(M2_2_PHEN_fit$estimate)
  AIC_results[6] = -2 * M3_1_fit$maximum + 2 * length(M3_1_fit$estimate)
  AIC_results[7] = -2 * M3_2_fit$maximum + 2 * length(M3_2_fit$estimate)
  AIC_results[8] = -2 * M3_3_fit$maximum + 2 * length(M3_3_fit$estimate)
  AIC_results[9] = -2 * M3_4_fit$maximum + 2 * length(M3_4_fit$estimate)
  AIC_results[10] = -2 * M3_5_fit$maximum + 2 * length(M3_5_fit$estimate)
  AIC_results[11] = -2 * M3_6_fit$maximum + 2 * length(M3_6_fit$estimate)
  
  
  AIC_results[12] = -2 * M4_fit$maximum + 2 * length(M4_fit$estimate)
  AIC_results[13] = -2 * M5_1_fit$maximum + 2 * length(M5_1_fit$estimate) 
  AIC_results[14] = -2 * (M5_2_SA_PA_fit$maximum + M5_2_SA_PB_fit$maximum + 
                           M5_2_SB_PA_fit$maximum + M5_2_SB_PB_fit$maximum) + 
    2 * 4 * length(M5_2_SA_PA_fit$estimate)
  
  AIC_results[15] = -2 * M5_2_PRED_fit$maximum + 2 * length(M5_2_PRED_fit$estimate)
  AIC_results[16] = -2 * M5_2_PHEN_fit$maximum + 2 * length(M5_2_PHEN_fit$estimate)
  AIC_results[17] = -2 * M6_1_fit$maximum + 2 * length(M6_1_fit$estimate)
  AIC_results[18] = -2 * M6_2_fit$maximum + 2 * length(M6_2_fit$estimate)
  AIC_results[19] = -2 * M6_3_fit$maximum + 2 * length(M6_3_fit$estimate)
  AIC_results[20] = -2 * M6_4_fit$maximum + 2 * length(M6_4_fit$estimate)
  AIC_results[21] = -2 * M6_5_fit$maximum + 2 * length(M6_5_fit$estimate)
  AIC_results[22] = -2 * M6_6_fit$maximum + 2 * length(M6_6_fit$estimate)
  
  
  which((AIC_results - min(AIC_results) <= 6))
  # Perform the required QAIC calculations: 
  
  QAIC_results = rep(0, times = 6)
  
  saturated_log_likelihood = factorial(num_perch) * (truth_set ^ truth_set) * ((num_perch - truth_set) ^ (num_perch - truth_set)) / 
    (factorial(truth_set) * (factorial(num_perch - truth_set)) * (num_perch ^ num_perch))
  
  saturated_log_likelihood = sum(log(saturated_log_likelihood))
  
  df = length(x) - 2
  v_tilda = (2 / df) * (saturated_log_likelihood - M2_fit$maximum)
  
  QAIC_results[1] = - (2 / v_tilda) * M1_fit$maximum + 2 * length(M1_fit$estimate)
  QAIC_results[2] = - (2 / v_tilda) * M2_fit$maximum + 2 * length(M2_fit$estimate)
  QAIC_results[3] = - (2 / v_tilda) * M3_fit$maximum + 2 * length(M3_fit$estimate)
  QAIC_results[4] = - (2 / v_tilda) * M4_fit$maximum + 2 * length(M4_fit$estimate)
  QAIC_results[5] = - (2 / v_tilda) * M5_fit$maximum + 2 * length(M5_fit$estimate)
  QAIC_results[6] = - (2 / v_tilda) * M6_fit$maximum + 2 * length(M6_fit$estimate)
  
  
  
  # Instantiate a vector to keep track of the KLD for the i'th iteration:
  zth_results = rep(0, time = 6)
  c = 0 # Initialise the constant, c. 
  
  # Calculate a KLD value for the fitted models above: 
  for(z in 1:Z)
  {
    
    alpha_truth = p2(x) / phi
    beta_truth = (1 - p2(x)) / phi
    
    # Generate another set of "truth" data points (as per previous) and their 
    # associated probabilities. 
    
    x_valid = extraDistr::rbbinom(n = length(x), size = num_perch, 
                                  alpha = alpha_truth, beta = beta_truth)
    
    # Calculate the probabilities associated with the x values we have generated above. 
    # @note: extraDistr's dbbinom takes alpha and beta as vectors and gives no 
    # unexpected behaviour, provided the vectors are the same length. 
    # Consider this to be in development until results have been validated. 
    
    # Produce a set of "truth" data points:
    x_prob = extraDistr::dbbinom(x = x_valid, size = num_perch, alpha = alpha_truth, beta = beta_truth)
    
    #c = sum((x_prob) * log(x_prob)) # Expected log probability of the i'th outcome. The c term is common to all models. 
    c = c + sum(log(x_prob))
    
    
    # Calculate the KLD values themselves: 
    zth_results[1] = zth_results[1] + sum(log(x_prob)) - sum(log(dbinom(x = x_valid, size = num_perch, prob = exp(-tau * M1_fit$estimate[1]))))
    
    
    zth_results[2] = zth_results[2] + sum(log(x_prob)) - sum(log(dbinom(x = x_valid, size = num_perch, prob = exp((-tau * M2_fit$estimate[1]) /
                                                                                                                    (1 + M2_fit$estimate[2] * x))))) 
    
    
    zth_results[3] = zth_results[3] + sum(log(x_prob)) - sum(log(dbinom(x = x_valid, size = num_perch,
                                                                        prob = exp(-tau * M3_fit$estimate[1] + M3_fit$estimate[2] * x) /
                                                                          (1 + exp(-tau * M3_fit$estimate[1] + M3_fit$estimate[2] * x))))) 
    
    
    M4_p_bar = exp(-tau * M4_fit$estimate[1])
    M4_alpha = M4_p_bar /  M4_fit$estimate[2]
    M4_beta = (1 - M4_p_bar) / M4_fit$estimate[2]
    
    zth_results[4] = zth_results[4] + sum(log(x_prob)) - sum(log(extraDistr::dbbinom(x = x_valid, size = num_perch,
                                                                                     alpha = M4_alpha,
                                                                                     beta = M4_beta)))
    
    
    M5_p_bar = exp((-tau * M5_fit$estimate[1]) / (1 + M5_fit$estimate[2] * x))
    M5_alpha = M5_p_bar /  M5_fit$estimate[3]
    M5_beta = (1 - M5_p_bar) / M5_fit$estimate[3]
    
    zth_results[5] = zth_results[5] + sum(log(x_prob)) - sum(log(extraDistr::dbbinom(x = x_valid, size = num_perch,
                                                                                     alpha = M5_alpha,
                                                                                     beta = M5_beta))) 
    
    
    M6_p_bar = exp(-tau * M6_fit$estimate[1] + M6_fit$estimate[2] * x) / 
      (1 + exp(-tau * M6_fit$estimate[1] + M6_fit$estimate[2] * x))
    M6_alpha = M6_p_bar /  M6_fit$estimate[3]
    M6_beta = (1 - M6_p_bar) / M6_fit$estimate[3]
    
    zth_results[6] = zth_results[6] + sum(log(x_prob)) - sum(log(extraDistr::dbbinom(x = x_valid, size = num_perch,
                                                                                     alpha = M6_alpha,
                                                                                     beta = M6_beta)))
    
    
  }
  
  return_list = vector("list", length = 4) # Initialise a list to return from the function.
  
  # Save the return values in a set order:
  # @pos1: the mean KLD for the current fit.
  # @pos2: the AIC estimate for the current fit, calculated using the mean KLD values and estimate of the contstant, c.
  # @pos3: the AIC values for the current fit, calculated using model max log likelihood.
  
  return_list[[1]] = zth_results / Z
  return_list[[2]] = 2 * (zth_results / Z - c / Z)
  return_list[[3]] = AIC_results
  return_list[[4]] = QAIC_results
  
  return(return_list)
  
  
}


# Set up the parallel environment:
num_cores = parallel::detectCores()
cluster = parallel::makeCluster(num_cores, type = "PSOCK")
registerDoSNOW(cluster)

pb = txtProgressBar(max = num_reps, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)

parallel_results = foreach(i=1:num_reps, 
                           .options.snow = opts) %dopar% process_rep(i)

close(pb)
stopCluster(cluster)

KLD_parallel_results = purrr::map(parallel_results, 1)
AIC_estimate_parallel_results = purrr::map(parallel_results, 2)
AIC_true_parallel_results = purrr::map(parallel_results, 3)
QAIC_true_parallel_results = purrr::map(parallel_results, 4)


KLD_mat = matrix(unlist(KLD_parallel_results), ncol = 6, byrow = TRUE)
AIC_estimate_mat = matrix(unlist(AIC_estimate_parallel_results), ncol = 6, byrow = TRUE)
AIC_mat = matrix(unlist(AIC_true_parallel_results), ncol = 6, byrow = TRUE)
QAIC_mat = matrix(unlist(QAIC_true_parallel_results), ncol = 6, byrow = TRUE)

# xlsx::write.xlsx(AIC_mat, file = "binomial_scenario_c.xlsx", sheetName = "AIC", row.names = FALSE, col.names = FALSE)
# xlsx::write.xlsx(KLD_mat, file = "binomial_scenario_c.xlsx", sheetName = "KLD_estimates", row.names = FALSE, col.names = FALSE, append = TRUE)
# xlsx::write.xlsx(AIC_estimate_mat, file = "binomial_scenario_c.xlsx", sheetName = "AIC_KLD_estimates", row.names = FALSE, col.names = FALSE, append = TRUE)
# xlsx::write.xlsx(QAIC_mat, file = "binomial_scenario_c.xlsx", sheetName = "QAIC", row.names = FALSE, col.names = FALSE, append = TRUE)

