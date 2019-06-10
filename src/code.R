rm(list = ls())
gc(reset = TRUE)

if(!require(tidyverse)) install.packages('tidyverse'); library(tidyverse)

group_index = c(1, 1, 2, 1, 2, 1, 2, 3, 3, 3)

set.seed(1)

simul_mat = matrix(rnorm(1000, 0.01, 0.1), nrow = 100)

#### hyperparameters

tau = 0.5; rho = 0.5

lambda_1 = 0.5; lambda_2 = 1

#### unchangeable

mu = colMeans(simul_mat)

mu_zero = 0.01

n = nrow(simul_mat); p = ncol(simul_mat)

Amat = rbind(diag(p), diag(p), 1, mu)

Bmat = rbind(-diag(2*p), 0, 0)

Cmat = c(rep(0, 2*p), 1, mu_zero)

X_tilde = t(cbind(-1, simul_mat))

K = cbind(0, diag(p)) 

pl = group_index %>% table %>% sqrt %>% as.vector()

lr1 = 0.01

lr2 = 0.01

#### initial value

tmp_z = rnorm(2*p, 0, 0.1)

tmp_u = rnorm(2*p + 2, 0, 0.01)

tmp_beta = rep(0, p)

tmp_beta_0 = 0

tmp_beta_tilde = c(tmp_beta_0, tmp_beta)

for(j in seq_len(10000)){
  
  #### step 1
  i = 1 
  
  while(i <= 10){  
    tmp_W = if_else(abs(((t(X_tilde) %*% tmp_beta_tilde)*4)) <= 1e-5, 100000, abs(1/((t(X_tilde) %*% tmp_beta_tilde)*4))) %>% as.vector() %>% diag()
    
    tmp_beta_tilde = (-1/2) * solve(X_tilde %*% tmp_W %*% t(X_tilde) + (rho/2) * t(Amat %*% K) %*% (Amat %*% K))  %*%
      ((tau - 1/2) * X_tilde %*% rep(1, n) + rho * t(Amat %*% K) %*% (Bmat %*% tmp_z) - rho * t(Amat %*% K) %*% Cmat + t((Amat %*% K)) %*% tmp_u )
    
    i = i + 1
  }
  #### step 2
  
  tmp_beta = tmp_beta_tilde[-1]
  
  v = (Amat %*% tmp_beta - Cmat)
  
  v1 = v[1:p]
  
  tmp_z1 = if_else(abs(tmp_z[1:p] - lr1 * (rho * (tmp_z[1:p] - v1) - tmp_u[1:p])) < lambda_1,  0,
                   tmp_z[1:p] - lr1 * (rho * (tmp_z[1:p] - v1) - tmp_u[1:p]) - lambda_1 * sign(tmp_z[1:p] - lr1 * (rho * (tmp_z[1:p] - v1) - tmp_u[1:p]) - lambda_1))
  
  #### step 3
  
  v2 = v[(p+1):(2*p)]
  u2 = tmp_u[(p+1):(2*p)]
  
  tmp_z2 = tmp_z[(p+1):(2*p)] - lr2 * (rho * (tmp_z[(p+1):(2*p)] - v2) - u2)
  
  for (l in unique(group_index)){
    pll = pl[l]
    v2l = v2[group_index == l]
    u2l = u2[group_index == l]
    
    if(sqrt(t(tmp_z2[group_index == l])%*%(tmp_z2[group_index == l])) > lambda_2 * sqrt(pll)){
      tmp_z2[group_index == l] = tmp_z2[group_index == l] - lambda_2 * sqrt(pll) * tmp_z2[group_index == l] / c(sqrt(t(tmp_z2[group_index == l])%*%(tmp_z2[group_index == l])))
    } else {tmp_z2[group_index == l] = 0 }
  }
  
  tmp_z = c(tmp_z1, tmp_z2)
  
  #### step 4
  
  tmp_u = tmp_u + rho * (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)
  
  loss = sum(tau * (t(X_tilde) %*% tmp_beta_tilde) * I(t(X_tilde) %*% tmp_beta_tilde >= 0) + (1-tau)*(-t(X_tilde) %*% tmp_beta_tilde)*I(t(X_tilde) %*% tmp_beta_tilde < 0)) + 
    lambda_1 * sum(abs(tmp_z1)) + 
    lambda_2 * sum(pl * tapply(tmp_z2, group_index, function(x) sqrt(t(x) %*% (x)))) + 
    t(tmp_u) %*% (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) + 
    (rho/2) * t(Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) %*% (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)
  
  #### Loss function
  if( j %% 100 == 0){
    cat(j, 'loss:: ', loss, 'dual term:: ',  t(tmp_u) %*% (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat),'\n')
    }
}

tmp_beta
tmp_z1
tmp_z2
tmp_u * (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) 
mu %*% tmp_beta

# primal feasibility 
Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat

sum(tmp_beta)

