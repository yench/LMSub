####################################################################################
# Purpose: functions for influence function-based variance estimation and functions
#          adapted from Li et al. (2023) for generating data satisfying landmark
#          models.
# Author: Yen Chang
# Creation Date: 3/22/25
####################################################################################

# functions to compute joint inclusion probabilities
library(Rcpp)
library(RcppArmadillo)
sourceCpp('Rho.cpp')

jointVP = function(datph2, m, wgtcol){
  ctrl.n = sum(datph2$status==0)
  case.times = datph2[status==1, T.obs]
  ctrl.times = datph2[status==0, T.obs]
  case.Z3 = datph2[status==1, Z3]
  ctrl.Z3 = datph2[status==0, Z3]
  
  H = datph2[status==1, (1-2*..m/(nriskZ3-1)+..m*(..m-1)/((nriskZ3-1)*(nriskZ3-2)))/(1-..m/(nriskZ3-1))^2]
  
  Rho = Rho(H, ctrl.n, case.times, ctrl.times, case.Z3, ctrl.Z3)
  Rhony = Rho + t(Rho)
  p0 = 1/datph2[status==0, ][[wgtcol]]
  
  Vij = Rhony * tcrossprod(1-p0) + diag((1-p0)*(p0))
  Pij = Rhony * tcrossprod(1-p0) + tcrossprod(p0)
  diag(Pij) = p0
  
  return(Vij/Pij)
}

# function to compute IF-based SE following Shin et al. (2020)
sqdg = function(x){sqrt(diag(x))}

# function to estimate increments in cumulative baseline hazards and their influence functions
sourceCpp('dLambda0.cpp')

#======================================================================
# Functions adapted from the R code accompanying Li et al. (2023) 
# (https://github.com/liwh0904/Compare_JM_and_LM)
#======================================================================
# PDF of theta(tj)
beta.distribution =  function(theta.var, j){
  gamma(a_beta_theta[j] + b_beta_theta[j])/(gamma(a_beta_theta[j]) * 
  gamma(b_beta_theta[j])) * theta.var^(a_beta_theta[j]-1) * (1-theta.var)^(b_beta_theta[j]-1)
}

# pre-defined Weibull baseline H_0(u, 0)
H_0_base <- function(t) {(lambda*t)^kappa}

# function to get H_0(u, tj) 
library(nleqslv)
H_0_solve = function(u, j){
  H_0_tj_u = sum(exp(-H_0_base(landmark.time[j] + u) * theta) * theta.pdf.input[1,] * delta.theta)
  H_0_tj = sum(exp(-H_0_base(landmark.time[j]) * theta) * theta.pdf.input[1,] * delta.theta)
  H_0_general = function(H_0){
    sum(exp(-H_0 * theta) * theta.pdf.input[j,] * delta.theta) - H_0_tj_u/H_0_tj
  }
  
  H_0_start = 1  ## initial value
  return(nleqslv(H_0_start, H_0_general)$x)
}