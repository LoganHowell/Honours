library(tidyverse)
library(Rfast)
library(epitools)
library(extraDistr)
library(foreach)
library(doParallel)
library(parallel)
library(stats4)

# Binomial Example: 

r_pool_replicates = 5

kelp_bunches = seq(1:4) # The number of kelp bunches. 
tau = 15 # The duration of the experiment. 
alpha = 0.075 # Pre-defined capture rate (per hour). 
beta = 0.5 # The effectiveness of the kelp in mitigating perch capture. 
num_truths = 10000 # The number of truth values to be generated. 
num_perch = 20 # The number of prey fish in the tank. 
phi = 0.1 # The pre-determined overdispersion factor. 


# Function for the true capture rate. 
# @param x: The number of bunches in the pool. 
true_capture_rate = function(x) {
  return(alpha / (1 + beta * x))
}

# Collection of functions associated with each of the mean probabilities of survival 
# for each model.
# @param x: The number of bunches in the pool (in each relevant function, below).
p1 = function() { 
  return(exp(-tau * alpha))  
}

p2 = function(x) {
  return(exp(-tau * alpha / (1 + beta * x)))
}

p3 = function(x) {
  return(exp((-tau * alpha + beta *x) / (1 + exp(-tau * alpha + beta * x))))
}

func_list = list(p1, p2, p3)


# Generate a list of x-values corresponding to the different treatments: 
x = rep(1:4, each = num_truths)

# Generate the "truth" dataset: 
truth_set = c()

for(i in 1:length(x))
{
  alpha_truth = p2(x[i]) / phi
  beta_truth = (1 - p2(x[i])) / phi
  
  truth_set[i] = extraDistr::rbbinom(n = 1, size = num_perch, 
                                     alpha = alpha_truth, beta = beta_truth)
}


test = function(x) { 
  return(-(x-2)^2)  
}


log_lik_theta_given_y = function(function_index, n, y) {
  
  # Get the relevant version of p_bar:
  if(function_index == 1) {
    p_bar = func_list[[function_index]]()
  } else
  {
    p_bar = func_list[[function_index]](y)
  }
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  return(Lgamma(num_perch + 1) + Lgamma(a + b) + Lgamma(y + a) + Lgamma(n - y + b) - 
    Lgamma(y +  1) - Lgamma(n - y + 1) - Lgamma(a) - Lgamma(b) - Lgamma(n + a + b))
}


LL_M4 = function(param, y, n) {
  
  p_bar = exp(-tau * param[1])
  
  a = p_bar / param[2]
  b = (1 - p_bar) / param[2]
  
  return(Lgamma(num_perch + 1) + Lgamma(a + b) + Lgamma(y + a) + Lgamma(n - y + b) - 
           Lgamma(y +  1) - Lgamma(n - y + 1) - Lgamma(a) - Lgamma(b) - Lgamma(n + a + b))
}


n = 20
y = 3
LL_M4 = function(alpha, phi, n, y) {
  
  p_bar = exp(-tau * alpha)
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  # print(paste("a:", a, "b:", b, sep = " "))
  # print(paste("alpha:", alpha, "phi:", phi, sep = "   "))
  # 
  # print(paste("Result:", -(Lgamma(n + 1) + Lgamma(a + b) + Lgamma(y + a) + Lgamma(n - y + b) - 
  #               Lgamma(y +  1) - Lgamma(n - y + 1) - Lgamma(a) - Lgamma(b) - Lgamma(n + a + b))))
  
  return(-(Lgamma(n + 1) + Lgamma(a + b) + Lgamma(y + a) + Lgamma(n - y + b) - 
           Lgamma(y +  1) - Lgamma(n - y + 1) - Lgamma(a) - Lgamma(b) - Lgamma(n + a + b)))
}

stats4::mle((LL_M4), start = list(alpha = 1, phi = 1), method = "L-BFGS-B", 
            lower = c(1e-10, 1e-10), upper = c(Inf, Inf), ...n = 3, y = 20)
