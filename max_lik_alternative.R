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

y = truth_set
num_perch = 20
tau = 15

LL_M2 = function(pars) {
  
  # print(paste("alpha: ", alpha))
  # print(paste("beta: ", beta))
  # print(paste("Result: ", (((-y * tau * alpha) / (1 + beta * y)) + log(1 - exp((-tau * alpha) / (1 + beta * y))) * (num_perch - y))))
  # print("")
  
  return(-(((-y * tau * pars[1]) / (1 + pars[2] * y)) + log(1 - exp((-tau * pars[1]) / (1 + pars[2] * y))) * (num_perch - y)))
}

maxLik::maxLik(logLik = LL_M2, start = c(1,1))






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


A_LL_M1 = matrix(c(1), 1, 1)
B_LL_M1 = matrix(c(0), 1, 1)

A_LL_M2 = matrix(c(1,0,0,1), 2,2)
B_LL_M2 = matrix(c(0,0), 2, 1)

A_LL_M4 = matrix(c(1,0,0,1), 2,2)
B_LL_M4 = matrix(c(0,0), 2, 1)

A_LL_M5 = matrix(c(1,0,0,0,1,0,0,0,1), 3,3)
B_LL_M5 = matrix(c(0,0, 0), 3, 1)

M1_fit = maxLik::maxLik(logLik = LL_M1, start = c(1), constraints = list(ineqA = A_LL_M1, ineqB = B_LL_M1), method = "BFGS")
M2_fit = maxLik::maxLik(logLik = LL_M2, start = c(1,1), constraints = list(ineqA = A_LL_M2, ineqB = B_LL_M2))
M3_fit = maxLik::maxLik(logLik = LL_M3, start = c(1,1))
M4_fit = maxLik::maxLik(logLik = LL_M4, start = c(1,1), constraints = list(ineqA = A_LL_M4, ineqB = B_LL_M4))
M5_fit = maxLik::maxLik(logLik = LL_M5, start = c(1,1,1), constraints = list(ineqA = A_LL_M5, ineqB = B_LL_M5))
M6_fit = maxLik::maxLik(logLik = LL_M6, start = c(1,1,1))




