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

# Binomial Example: 

r_pool_replicates = 5
kelp_bunches = seq(1:4) # The number of kelp bunches. 
tau = 15 # The duration of the experiment. 
alpha = 0.075 # Pre-defined capture rate (per hour). 
beta = 0.5 # The effectiveness of the kelp in mitigating perch capture. 
num_truths = 5 # The number of truth values to be generated and which the models will be fit to (for each treatment). 
num_perch = 20 # The number of prey fish in the tank. 
phi = 0.1 # The pre-determined overdispersion factor.
num_reps = 100 # The number of times the KLD will be calculated. 
Z = 10000 # The number of new points to use in calculation of cumulative IC (for each treatment). 


# Constraint sets for later optimisation:
A_LL_M1 = matrix(c(1), 1, 1)
B_LL_M1 = matrix(c(0), 1, 1)

A_LL_M2 = matrix(c(1,0,0,1), 2,2)
B_LL_M2 = matrix(c(0,0), 2, 1)

A_LL_M4 = matrix(c(1,0,0,1), 2,2)
B_LL_M4 = matrix(c(0,0), 2, 1)

A_LL_M5 = matrix(c(1,0,0,0,1,0,0,0,1), 3,3)
B_LL_M5 = matrix(c(0,0, 0), 3, 1)


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


# Generate a list of x-values corresponding to the different treatments: 
x = rep(1:4, each = num_truths)


# We remove the log of the combination choose(num_perch, y) within each of M1 - M3
# as it is constant with respect to alpha in the maximisation of the log likelihood. 
LL_M1 = function(pars) {
  
  return((-tau * pars[1] * y + (num_perch - y) * log(1 - exp(-tau * pars[1]))))
}

LL_M2 = function(pars) {
  
  # print(paste("pars[1]: ", pars[1]))
  # print(paste("pars[2]: ", pars[2]))
  # print(paste("Result: ", (((-y * tau * pars[1]) / (1 + pars[2] * y)) + log(1 - exp((-tau * pars[1]) / (1 + pars[2] * y))) * (num_perch - y))))
  # print("")
  
  return(sum(((-y * tau * pars[1]) / (1 + pars[2] * y)) + log(1 - exp((-tau * pars[1]) / (1 + pars[2] * y))) * (num_perch - y)))
}

LL_M3 = function(par) {
  
  numerator = exp(-tau * par[1] + par[2] * y)
  denominator = 1 + exp(-par[1] * tau + par[2] * y)
  
  # print(paste("Numerator::", numerator, sep = "  "))
  # print(paste("Denominator::", denominator, sep = "  "))
  # 
  # print((y * (log(numerator) - log(denominator)) + (num_perch - y) * log(1 - (numerator / denominator))))
  
  return(sum(y * (log(numerator) - log(denominator)) + (num_perch - y) * log(1 - (numerator / denominator))))
}


LL_M4 = function(par) {
  
  p_bar = exp(-tau * par[1])
  
  a = p_bar / par[2]
  b = (1 - p_bar) / par[2]
  
  # print(paste("a:", a, "b:", b, sep = " "))
  # print(paste("par[1]:", par[1], "par[2]:", par[2], sep = "   "))
  # 
  # print(paste("Result:", -sum((lgamma(n + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(n - y + b) -
  #               lgamma(y +  1) - lgamma(n - y + 1) - lgamma(a) - lgamma(b) - lgamma(n + a + b)))))
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(num_perch - y + b) - 
               lgamma(y +  1) - lgamma(num_perch - y + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


LL_M5 = function(par) {
  
  p_bar = exp(-tau * par[1] / (1 + par[2] * y))
  
  a = p_bar / par[3]
  b = (1 - p_bar) / par[3]
  
  # print(paste("a:", a, "b:", b, sep = " "))
  # print(paste("par[1]:", par[1], "par[3]:", par[3], sep = "   "))
  # 
  # print(paste("Result:", -(lgamma(n + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(n - y + b) - 
  #               lgamma(y +  1) - lgamma(n - y + 1) - lgamma(a) - lgamma(b) - lgamma(n + a + b))))
  
  return(sum(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(num_perch - y + b) - 
               lgamma(y +  1) - lgamma(num_perch - y + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
}


LL_M6 = function(par) {
  
  p_bar = exp(-tau * par[1] + par[2] * y) / (1 + exp(-tau * par[1] + par[2] * y))
  
  a = p_bar / par[3]
  b = (1 - p_bar) / par[3]
  
  # print(paste("p_bar:", p_bar))
  # print(paste("a:", a, "b:", b, sep = " "))
  # print(paste("par[1]:", par[1], "par[3]:", par[3], sep = " "))
  # print(paste("Result:", -(lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(num_perch - y + b) -
  #                            lgamma(y +  1) - lgamma(num_perch - y + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b))))
  # print("")
  
  
  return(sum((lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(num_perch - y + b) - 
                lgamma(y +  1) - lgamma(num_perch - y + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b))))
  
}

ll_mat = matrix(data = 0, nrow = num_reps, ncol = 6)

for(i in 1:num_reps)
{
  # Generate the i'th "truth" dataset: 
  truth_set = c()
  
  for(j in 1:length(x))
  {
    p_bar = p2(x[j])
    
    alpha_truth = p_bar / phi
    beta_truth = (1 - p_bar) / phi
    
    truth_set[j] = extraDistr::rbbinom(n = 1, size = num_perch, 
                                       alpha = alpha_truth, beta = beta_truth)
  }
  
  # Dear Lord, please forgive me for what I'm about to do:
  # Set y to be the current truth set all optimisation functions have global access :_( 
  y = truth_set 
  
  print(i)
  
  # Fit each of the models to the generated data:
   
  M1_fit = maxLik::maxLik(logLik = LL_M1, start = c(1), constraints = list(ineqA = A_LL_M1, ineqB = B_LL_M1), method = "BFGS")
  M2_fit = maxLik::maxLik(logLik = LL_M2, start = c(1,1), constraints = list(ineqA = A_LL_M2, ineqB = B_LL_M2))
  M3_fit = maxLik::maxLik(logLik = LL_M3, start = c(1,1))
  M4_fit = maxLik::maxLik(logLik = LL_M4, start = c(1,1), constraints = list(ineqA = A_LL_M4, ineqB = B_LL_M4))
  M5_fit = maxLik::maxLik(logLik = LL_M5, start = c(1,1,1), constraints = list(ineqA = A_LL_M5, ineqB = B_LL_M5))
  M6_fit = maxLik::maxLik(logLik = LL_M6, start = c(1,1,1))
  
  # Instantiate a vector to keep track of the KLD for the i'th iteration:
  zth_results = rep(0, time = 6)
  
  # Calculate a KLD value for the fitted models above: 
  for(z in 1:Z)
  {
    
    # The vector of validation observations and their associated probabilities: 
    x_valid = rep(0, times = length(x))

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
    x_prob = extraDistr::dbbinom(x = x_valid, size = num_perch, alpha = alpha_truth, beta = beta_truth)
    
    # Calculate the KLD values themselves: 
    zth_results[1] = zth_results[1] + (sum(log(dbinom(x = x_valid, size = num_perch, prob = exp(-tau * M1_fit$estimate[1])))) - 
      sum(log(x_prob))) / Z
    
    zth_results[2] = zth_results[2] + (sum(log(dbinom(x = x_valid, size = num_perch, prob = exp((-tau * M2_fit$estimate[1]) /
                                                                      (1 + M2_fit$estimate[2] * x_valid))))) - 
      sum(log(x_prob))) / Z
    
    zth_results[3] = zth_results[3] + (sum(log(dbinom(x = x_valid, size = num_perch,
               prob = exp(-tau * M3_fit$estimate[1] + M3_fit$estimate[2] * x_valid) /
                 (1 + exp(-tau * M3_fit$estimate[1] + M3_fit$estimate[2] * x_valid))))) - 
      sum(log(x_prob))) / Z
    
    zth_results[4] = zth_results[4] + (sum(log(extraDistr::dbbinom(x = x_valid, size = num_perch,
                            alpha = M4_fit$estimate[1],
                            beta = M4_fit$estimate[2]))) - 
      sum(log(x_prob))) / Z
    
    zth_results[5] = zth_results[5] + (sum(log(extraDistr::dbbinom(x = x_valid, size = num_perch,
                            alpha = M5_fit$estimate[1],
                            beta = M5_fit$estimate[2]))) - 
      sum(log(x_prob))) / Z
    
    zth_results[6] = zth_results[6] + (sum(log(extraDistr::dbbinom(x = x_valid, size = num_perch,
                            alpha = M6_fit$estimate[1],
                            beta = M6_fit$estimate[2]))) - 
      sum(log(x_prob))) / Z
    
  }
  
  ll_mat[i, ] = zth_results

}

EKLD = colMeans(ll_mat)


