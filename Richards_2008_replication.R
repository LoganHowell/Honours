library(tidyverse)
library(Rfast)
library(epitools)
library(extraDistr)
library(foreach)
library(doParallel)
library(parallel)
library(stats4)
library(fAsianOptions)

# Binomial Example: 

r_pool_replicates = 5

kelp_bunches = seq(1:4) # The number of kelp bunches. 
tau = 15 # The duration of the experiment. 
alpha = 0.075 # Pre-defined capture rate (per hour). 
beta = 0.5 # The effectiveness of the kelp in mitigating perch capture. 
num_truths = 500 # The number of truth values to be generated and which the models will be fit to (for each treatment). 
num_perch = 20 # The number of prey fish in the tank. 
phi = 0.1 # The pre-determined overdispersion factor. 
Z = 500 # The number of new points to use in calculation of cumulative IC (for each treatment). 


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


# We remove the log of the combination choose(num_perch, y) within each of M1 - M3
# as it is constant with respect to alpha in the maximisation of the log likelihood. 
LL_M1 = function(alpha) {
  
  return(-(-tau * alpha * y + (num_perch - y) * log(1 - exp(-tau * alpha))))
}


LL_M2 = function(alpha, beta) {
  
  return(-(((-y * tau * alpha) / (1 + beta * y)) + log(1 - exp((-tau * alpha) / (1 + beta * y))) * (num_perch - y)))
}


LL_M3 = function(alpha, beta) {
  
  numerator = exp(-tau * alpha + beta * y)
  denominator = 1 + exp(-alpha * tau + beta * y)
  
  # print(paste("Numerator::", numerator, sep = "  "))
  # print(paste("Denominator::", denominator, sep = "  "))
  # 
  # print((y * (log(numerator) - log(denominator)) + (num_perch - y) * log(1 - (numerator / denominator))))
  
  return(-(y * (log(numerator) - log(denominator)) + (num_perch - y) * log(1 - (numerator / denominator))))
}

LL_M4 = function(alpha, phi) {
  
  p_bar = exp(-tau * alpha)
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  # print(paste("a:", a, "b:", b, sep = " "))
  # print(paste("alpha:", alpha, "phi:", phi, sep = "   "))
  # 
  # print(paste("Result:", -(lgamma(n + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(n - y + b) - 
  #               lgamma(y +  1) - lgamma(n - y + 1) - lgamma(a) - lgamma(b) - lgamma(n + a + b))))
  
  return(-(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(num_perch - y + b) - 
             lgamma(y +  1) - lgamma(num_perch - y + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


LL_M5 = function(alpha, beta, phi) {
  
  p_bar = exp(-tau * alpha / (1 + beta * y))
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  # print(paste("a:", a, "b:", b, sep = " "))
  # print(paste("alpha:", alpha, "phi:", phi, sep = "   "))
  # 
  # print(paste("Result:", -(lgamma(n + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(n - y + b) - 
  #               lgamma(y +  1) - lgamma(n - y + 1) - lgamma(a) - lgamma(b) - lgamma(n + a + b))))
  
  return(-(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(num_perch - y + b) - 
             lgamma(y +  1) - lgamma(num_perch - y + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}

LL_M6 = function(alpha, beta, phi) {
  
  p_bar = exp(-tau * alpha + beta * y) / (1 + exp(-tau * alpha + beta * y))
  
  a = p_bar / phi
  b = (1 - p_bar) / phi
  
  # print(paste("p_bar:", p_bar))
  # print(paste("a:", a, "b:", b, sep = " "))
  # print(paste("alpha:", alpha, "phi:", phi, sep = " "))
  # print(paste("Result:", -(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(num_perch - y + b) -
  #                            lgamma(y +  1) - lgamma(num_perch - y + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b))))
  # print("")
  
  
  return(-(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(num_perch - y + b) - 
             lgamma(y +  1) - lgamma(num_perch - y + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


# Access method for optimised parameters:
#opt_obj@details$par

for(i in 1:length(truth_set))
{
  # Dear Lord, please forgive me for what I'm about to do:
  # Extract the i'th value and overwrite the previously declared global variable 
  # so all optimisation functions have global access :_( 
  y = truth_set[i] 
  
  print(i)
  
  # Fit each of the models using MLE:
  if(y != 20 & y != 0)
  {
    M1_fit = stats4::mle((LL_M1), start = list(alpha = 1), method = "L-BFGS-B", 
                         lower = c(1e-3), upper = c(Inf))
  }
  
  # M2_fit = stats4::mle((LL_M2), start = list(alpha = 1, beta = 1), method = "L-BFGS-B", 
  #                      lower = c(1e-12, 1e-12), upper = c(Inf, Inf))
  # 
  # if(y != 2)
  # {
  #   M3_fit = stats4::mle((LL_M3), start = list(alpha = 1, beta = 1), method = "L-BFGS-B", 
  #                        lower = c(-Inf, -Inf), upper = c(Inf, Inf))
  # }
  # 
  # 
  # if(y != 9 & y != 5)
  # {
  #   M4_fit = stats4::mle((LL_M4), start = list(alpha = 1, phi = 1), method = "L-BFGS-B", 
  #                        lower = c(1e-12, 1e-12), upper = c(Inf, Inf))
  # }
  # 
  # M5_fit = stats4::mle((LL_M5), start = list(alpha = 1, beta =1, phi = 1), method = "L-BFGS-B", 
  #             lower = c(1e-12, 1e-12, 1e-12), upper = c(Inf, Inf, Inf))
  # 
  # if(y != 5 & y != 13 & y != 16)
  # {
  #   M6_fit = stats4::mle((LL_M6), start = list(alpha = 1, beta = 1, phi = 1), method = "L-BFGS-B", 
  #                        lower = c(-Inf, -Inf, 1e-12), upper = c(Inf, Inf, Inf))
  # }
  
  IC_M1 = rep(0, times = length(truth_set))
  IC_M2 = rep(0, times = length(truth_set))
  IC_M3 = rep(0, times = length(truth_set))
  IC_M4 = rep(0, times = length(truth_set))
  IC_M5 = rep(0, times = length(truth_set))
  IC_M6 = rep(0, times = length(truth_set))
  
  r_valid = rep(1:4, each = Z)
  
  # Generate the "truth" dataset: 
  truth_set_valid = c()
  true_probs = c()
  
  for(i in 1:length(r_valid))
  {
    alpha_truth = p2(r_valid[i]) / phi
    beta_truth = (1 - p2(r_valid[i])) / phi
    
    truth_set_valid[i] = extraDistr::rbbinom(n = 1, size = num_perch, 
                                       alpha = alpha_truth, beta = beta_truth)
    
    true_probs[i] = dbbinom(x = truth_set_valid[i], size = num_perch, alpha = alpha_truth, beta = beta_truth)
  }
  
  for(z in 1:length(truth_set_valid))
  {
    IC_M1[z] = IC_M1[z] + log(dbinom(x = truth_set_valid[z], size = num_perch, prob = exp(-tau * M1_fit@details$par[["alpha"]]))) - 
      log(true_probs[z])
    
    # IC_M2[z] = IC_M2[z] + log(dbinom(x = truth_set[z], size = num_perch, prob = exp((-tau * M2_fit@details$par[["alpha"]]) / 
    #                                                                                   (1 + M2_fit@details$par[["beta"]] * truth_set[z])))) -
    #   log(true_probs[z])
    # 
    # IC_M3[z] = IC_M3[z] + log(dbinom(x = truth_set[z], size = num_perch, 
    #                                  prob = exp(-tau * M3_fit@details$par[["alpha"]] + M3_fit@details$par[["beta"]] * truth_set[z]) /
    #                                    (1 + exp(-tau * M3_fit@details$par[["alpha"]] + M3_fit@details$par[["beta"]] * truth_set[z])))) - 
    #   log(true_probs[z])
      
  }
  
  IC_M1[i] = IC_M1[i] / length(truth_set_valid)
  IC_M2[i] = IC_M2[i] / length(truth_set_valid)
  IC_M3[i] = IC_M3[i] / length(truth_set_valid)
  
}

IC_M1_mean = mean(IC_M1)

