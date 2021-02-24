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

y = 10
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
  
  
  return((lgamma(num_perch + 1) + lgamma(a + b) + lgamma(y + a) + lgamma(num_perch - y + b) - 
             lgamma(y +  1) - lgamma(num_perch - y + 1) - lgamma(a) - lgamma(b) - lgamma(num_perch + a + b)))
  
}

LL_M1 = function(pars) {
  
  return((-tau * pars[1] * y + (num_perch - y) * log(1 - exp(-tau * pars[1]))))
}

y = 20
test = maxLik::maxLik(logLik = LL_M6, start = c(1,1,1))

A_LL_M1 = matrix(c(1), 1, 1)
B_LL_M1 = matrix(c(0), 1, 1)

test = maxLik::maxLik(logLik = LL_M1, start = c(1), constraints = list(ineqA = A_LL_M1, ineqB = B_LL_M1), method = "BFGS")
test$estimate
