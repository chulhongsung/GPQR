rm(list = ls())
gc(reset = TRUE)

if(!require(tidyverse)) install.packages('tidyverse'); library(tidyverse)

group_index = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)

set.seed(2)

simul_mat = matrix(c(rnorm(400, 0.1, 0.01), rnorm(100, 0.1, 0.01), rnorm(200, 0, 0.01), rnorm(300, 0, 0.01)), nrow = 100)

#### unchangeable

mu_zero = 0.1

mu = colMeans(simul_mat)

n = nrow(simul_mat); p = ncol(simul_mat)

Amat = rbind(diag(p), diag(p), 1, mu)

Bmat = rbind(-diag(2*p), 0, 0)

Cmat = c(rep(0, 2*p), 1, mu_zero)

X_tilde = t(cbind(-1, simul_mat))

K = cbind(0, diag(p)) 

pl = group_index %>% table %>% sqrt %>% as.vector()

tau = 0.5; rho = 0.7

lambda_1 = 1; lambda_2 = 20

#### initial value

tmp_z = rep(1, 2*p)/p

tmp_u = rep(0, 2*p+2)

tmp_beta = rep(1, p)/p

tmp_beta_0 = 0

tmp_beta_tilde = c(tmp_beta_0, tmp_beta)

norm_vec <- function(x) sqrt(sum(x^2))

loss_ = c()

dual_term = c()

for(j in seq_len(10000)){
  
  #### Loss function
  
  tmp_z1 = tmp_z[1:p]
  
  tmp_z2 = tmp_z[(p+1):(2*p)]
  
  loss = sum(tau * (t(X_tilde) %*% tmp_beta_tilde) * I(t(X_tilde) %*% tmp_beta_tilde >= 0) + (1-tau)*(-t(X_tilde) %*% tmp_beta_tilde)*I(t(X_tilde) %*% tmp_beta_tilde < 0)) + 
    lambda_1 * sum(abs(tmp_z1)) + 
    lambda_2 * sum(pl * tapply(tmp_z2, group_index, function(x) norm_vec(x))) + 
    t(tmp_u) %*% (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) + 
    (rho/2) * norm_vec(Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)^2
  
  loss_ = c(loss_, loss)
  
  dual_term = c(dual_term,  t(tmp_u) %*% (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat))
  
  if(j == 1){
    cat('initial loss::', loss, '\n')
    cat('==========================================================================================================', '\n')
  } 
  
  if(j %% 100 == 0){
    cat('iteration ', j, '\n')
    cat('loss:: ', loss, 'dual term:: ',  t(tmp_u) %*% (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat),'\n')
  }
  #### step 1
  i = 1
  
  while(i <= 10){
    
    tmp_W = if_else(abs(((t(X_tilde) %*% tmp_beta_tilde)*4)) <= 1e-5, 100000, abs(1/((t(X_tilde) %*% tmp_beta_tilde)*4))) %>% as.vector() %>% diag()
    
    tmp_beta_tilde = (-1/2) * solve(X_tilde %*% tmp_W %*% t(X_tilde) + (rho/2) * t(Amat %*% K) %*% (Amat %*% K)) %*%
      ((tau - 1/2) * X_tilde %*% rep(1, n) + rho * t(Amat %*% K) %*% ((1/rho) * tmp_u + Bmat %*% tmp_z - Cmat))
    
    i = i +1
  }
  
  #### step 2
  gradient_beta = tau * X_tilde %*% I(t(X_tilde) %*% tmp_beta_tilde >= 0) - (1 - tau) * X_tilde %*% I(t(X_tilde) %*% tmp_beta_tilde < 0) + t(Amat %*% K) %*% tmp_u +
    rho * t(Amat %*% K) %*% ((Amat %*% K) %*% tmp_beta_tilde + Bmat %*% tmp_z - Cmat)
  
  if(j %% 100 == 0){
    cat('KKT condition stationarity(beta)::  ', abs(gradient_beta) <= 1e-4 ,'\n')
  }
  
  tmp_beta = tmp_beta_tilde[-1]
  
  v = (Amat %*% tmp_beta - Cmat)
  
  v1 = v[1:p]
  
  tmp_z1 = v1 + tmp_u[1:p]/rho
  
  tmp_z1 = if_else(abs(tmp_z1) <= lambda_1, 0, tmp_z1 - lambda_1 * sign(tmp_z1 - lambda_1))
  
  if(j %% 100 == 0){
    cat('KKT condition stationarity(z1)::  ', abs(rho * (tmp_z1 - v1) - tmp_u[1:p]) <= lambda_1 ,'\n')
  }
  
  #### step 3
  
  v2 = v[(p+1):(2*p)]
  u2 = tmp_u[(p+1):(2*p)]
  
  tmp_z2 = v2 + u2/rho
  
  for (l in unique(group_index)){
    pll = pl[l]
    v2l = v2[group_index == l]
    u2l = u2[group_index == l]
    
    if(norm_vec(tmp_z2[group_index == l]) > lambda_2 * pll){
      tmp_z2[group_index == l] = tmp_z2[group_index == l] - lambda_2 * pll * tmp_z2[group_index == l] / norm_vec(tmp_z2[group_index == l])
    } else {tmp_z2[group_index == l] = 0 }
  }
  
  if(j %% 100 == 0){
    cat('KKT condition stationary(z2)::  ', tapply(tmp_z2, group_index, function(x) norm_vec(x)) <= lambda_2 ,'\n')
  }
  tmp_z = c(tmp_z1, tmp_z2)
  
  #### step 4
  tmp_u = tmp_u + rho * (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)
  
  if(j %% 100 == 0){
    cat('KKT condition dual feasibility::  ', abs(Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) <= 1e-4 ,'\n')
    cat('==========================================================================================================', '\n')
  }
}

tmp_beta

sum(tmp_beta)

mu %*% tmp_beta

tmp_z1

tmp_z2

