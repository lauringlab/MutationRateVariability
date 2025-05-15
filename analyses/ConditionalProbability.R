library(tidyverse)
library(extraDistr)

#helper functions to parameterize gamma-Poisson using
#mean and overdispersion index
shape_gpois <- function(mean, o_index){
  var <- o_index*mean
  return((mean^2)/(var-mean))
}

#probability of X = m given alive
p_m_given_alive <- function(m, m_mean, m_var){
  m_shape <- shape_gpois(m_mean, m_var/m_mean)
  dgpois(m, shape = m_shape, rate = m_shape/m_mean)
}

#probability of alive given X = m
p_alive_given_m <- function(m, p_lethal){
  return((1-p_lethal)^m)
}

#calculate the sum of P(M | alive) / P(alive | M)
#convergence condition is < 1e-7
calc_sum <- function(m_mean, m_var, p_lethal){
  p_m <- c(0,1)
  p_m <- p_m_given_alive(p_m, m_mean, m_var)/p_alive_given_m(p_m, p_lethal)
  i <- 3
  while(sum(p_m) - sum(p_m[-length(p_m)]) > 1e-7){
    p_m[i] <- p_m_given_alive(i-1, m_mean, m_var)/p_alive_given_m(i-1, p_lethal)
    #p_m[i] <- calc_num(i-1, m_mean, m_var, p_lethal)
    i <- i+1
  }
  return(sum(p_m))
}

#probability X = m, unconditional
p_m_uc <- function(m, m_mean, m_var, p_lethal){
  #P(M | alive) / P(alive | M)
  p_givens <- p_m_given_alive(m, m_mean, m_var)/p_alive_given_m(m, p_lethal)
  return(p_givens/calc_sum(m_mean, m_var, p_lethal))
}

#expected value of X, unconditional
exp_x_uc <- function(m_mean, m_var, p_lethal){
  if(m_var < m_mean){
    return(NA)
  }
  p_m <- c(0,1)
  p_m <- p_m * p_m_uc(p_m, m_mean, m_var, p_lethal)
  i <- 2
  while(sum(p_m) - sum(p_m[-length(p_m)]) > 1e-7){
    p_m[i] <- (i-1) * p_m_uc(i-1, m_mean, m_var, p_lethal)
    i <- i+1
  }
  return(sum(p_m))
}

#variance of X, unconditional
var_x_uc <- function(m_mean, m_var, p_lethal){
  if(m_var < m_mean){
    return(NA)
  }
  p_m <- c(0,1)
  p_m <- (p_m^2) * p_m_uc((p_m^1), m_mean, m_var, p_lethal)
  i <- 2
  while(sum(p_m) - sum(p_m[-length(p_m)]) > 1e-7){
    p_m[i] <- ((i-1)^2) * p_m_uc(i-1, m_mean, m_var, p_lethal)
    i <- i+1
  }
  return(sum(p_m) - (exp_x_uc(m_mean, m_var, p_lethal)^2))
}

#vectorize expected value and variance functions
exp_x_ucV <- Vectorize(exp_x_uc)
var_x_ucV <- Vectorize(var_x_uc)
