rm(list = ls())
gc(reset = TRUE)

if(!require(tidyverse)) install.packages('tidyverse'); library(tidyverse)

#### sgl.fit function

sgl.fit = function(data, group_index, mu_zero, tau, rho = 0.5, lambda_1, lambda_2, iter = 10000, verbose = TRUE){
  
  #### Design matrix & constraints
  
  mu = colMeans(data)
  n = nrow(data); p = ncol(data)
  Amat = rbind(diag(p), diag(p), 1, c(mu))
  Bmat = rbind(-diag(2*p), 0, 0)
  Cmat = c(rep(0, 2*p), 1, mu_zero)
  X_tilde = t(cbind(-1, data))
  K = cbind(0, diag(p)) 
  pl = group_index %>% table %>% sqrt %>% as.vector()
  
  #### Initial value
  
  tmp_z = rep(1, 2*p)/p
  tmp_u = rep(0, 2*p+2)
  tmp_beta = rep(1, p)/p
  tmp_beta_0 = 0
  tmp_beta_tilde = c(tmp_beta_0, tmp_beta)
  
  norm_vec <- function(x) sqrt(sum(x^2))
  
  j = 1
  
  while(j <= iter){
    
    #### Loss function
    
    tmp_z1 = tmp_z[1:p]
    tmp_z2 = tmp_z[(p+1):(2*p)]
    
    loss = sum(tau * (t(X_tilde) %*% tmp_beta_tilde) * I(t(X_tilde) %*% tmp_beta_tilde >= 0) + (1-tau)*(-t(X_tilde) %*% tmp_beta_tilde)*I(t(X_tilde) %*% tmp_beta_tilde < 0)) + 
      lambda_1 * sum(abs(tmp_z1)) + 
      lambda_2 * sum(pl * tapply(tmp_z2, group_index, function(x) norm_vec(x))) + 
      t(tmp_u) %*% (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) + 
      (rho/2) * norm_vec(Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)^2
    
    if(j == 1){
      cat('=========================================================', '\n')
      cat('initial loss: ', loss, '\n')
      cat('=========================================================', '\n')
    } 
    
    if((verbose != FALSE) & (j %% 1000 == 0)){
      cat('=========================================================', '\n')
      cat('iteration ', j, '\n')
      cat('loss: ', loss, '\n')
    }
    #### Step 1
    i = 1
    
    while(i <= 10){
      
      tmp_W = if_else(abs(((t(X_tilde) %*% tmp_beta_tilde)*4)) <= 1e-5, 100000, abs(1/((t(X_tilde) %*% tmp_beta_tilde)*4))) %>% as.vector() %>% diag()
      tmp_beta_tilde = (-1/2) * solve(X_tilde %*% tmp_W %*% t(X_tilde) + (rho/2) * t(Amat %*% K) %*% (Amat %*% K)) %*%
        ((tau - 1/2) * X_tilde %*% rep(1, n) + rho * t(Amat %*% K) %*% ((1/rho) * tmp_u + Bmat %*% tmp_z - Cmat))
      
      i = i +1
    }
    
    #### Step 2
    gradient_beta = tau * X_tilde %*% I(t(X_tilde) %*% tmp_beta_tilde >= 0) - (1 - tau) * X_tilde %*% I(t(X_tilde) %*% tmp_beta_tilde < 0) + t(Amat %*% K) %*% tmp_u +
      rho * t(Amat %*% K) %*% ((Amat %*% K) %*% tmp_beta_tilde + Bmat %*% tmp_z - Cmat)
    
    kkt1 = all(abs(gradient_beta) <= 1e-4)
    if((verbose != FALSE) & (j %% 1000 == 0)){
      cat('KKT condition stationarity1: ', kkt1, '\n')
    }
    
    tmp_beta = tmp_beta_tilde[-1]
    
    v = (Amat %*% tmp_beta - Cmat)
    v1 = v[1:p]
    
    tmp_z1 = v1 + tmp_u[1:p]/rho
    tmp_z1 = if_else(abs(tmp_z1) <= lambda_1, 0, tmp_z1 - lambda_1 * sign(tmp_z1 - lambda_1))
    
    kkt2 = all(abs(rho * (tmp_z1 - v1) - tmp_u[1:p]) <= lambda_1)
    
    if((verbose != FALSE) & (j %% 1000 == 0)){
      cat('KKT condition stationarity2: ', kkt2 ,'\n')
    }
    
    #### Step 3
    
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
    
    kkt3 = all(tapply(tmp_z2, group_index, function(x) norm_vec(x)) <= lambda_2)
    
    if((verbose != FALSE) & (j %% 1000 == 0)){
      cat('KKT condition stationarity3: ', kkt3, '\n')
    }
    
    tmp_z = c(tmp_z1, tmp_z2)
    
    #### Step 4
    tmp_u = tmp_u + rho * (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)
    
    kkt4 = all(abs(Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) <= 1e-4)
    
    if((verbose != FALSE) & (j %% 1000 == 0)){
      cat('KKT condition dual feasibility: ',  kkt4, '\n')
      cat('=========================================================', '\n')
    }
    
    if(all(c(kkt1, kkt2, kkt3, kkt4) == TRUE)){
      cat('=========================================================', '\n')
      cat('iteration ', j, 'loss: ', loss,' KKT condition satisfied!', '\n' )
      cat('KKT condition stationarity1: ', kkt1, '\n')
      cat('KKT condition stationarity2: ', kkt2 ,'\n')
      cat('KKT condition stationarity3: ', kkt3, '\n')
      cat('KKT condition dual feasibility: ',  kkt4, '\n')
      cat('=========================================================', '\n')
      
      return(list(solution = list(beta = c(tmp_beta_tilde), z = c(tmp_z), u = c(tmp_u)))); break
    }
    
    j =  j + 1
  }
  
  cat('Solution does not converges in', iter, 'iterations!', '\n')
  
  return(list(current_solution = list(beta = c(tmp_beta_tilde), z = c(tmp_z), u = c(tmp_u))))
}

#### Fit SGL

group_index = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3)

set.seed(1)

#### Simulation data matrix

simul_mat = matrix(c(rnorm(400, 0.1, 0.01), rnorm(100, 0.1, 0.01), rnorm(200, 0, 0.01), rnorm(300, 0, 0.01)), nrow = 100)

fit = sgl.fit(data=simul_mat, group_index=group_index, mu_zero=0.1, tau=0.5, rho=0.8, lambda_1=1, lambda_2=1)

#### Optimal solution
fit$solution
