if(!require(dplyr)) install.packages('dplyr'); library(dplyr)

#### sgl.fit function

pess_sgl_fit = function(data, group_index, mu_zero = NULL, tau, lambda_1, lambda_2, iter = 10000, verbose = TRUE, num = 1000){
  
  #### Design matrix & constraints
  norm_vec <- function(x) sqrt(sum(x^2))
  
  if(!is.null(mu_zero)){
    data = data * 100
    mu = colMeans(data)
    n = nrow(data); p = ncol(data)
    rho = 1
    Amat = rbind(diag(p), diag(p), 1, c(mu))
    Bmat = rbind(-diag(2*p), 0, 0)
    Cmat = c(rep(0, 2*p), 1, mu_zero)
    X_tilde = cbind(-1, data)
    K = cbind(0, diag(p)) 
    pl = group_index %>% table %>% sqrt %>% as.vector()
    
    #### Initial value
    
    tmp_z = rep(1, 2*p)/p
    tmp_u = rep(0, 2*p+2)
    tmp_beta = rep(1, p)/p
    tmp_beta_0 = 1e-5
    tmp_beta_tilde = c(tmp_beta_0, tmp_beta)
  }
  else {
    mu = colMeans(data)
    n = nrow(data); p = ncol(data)
    rho = 1
    Amat = rbind(diag(p), diag(p), 1)
    Bmat = rbind(-diag(2*p), 0)
    Cmat = c(rep(0, 2*p), 1)
    X_tilde = cbind(-1, data)
    K = cbind(0, diag(p)) 
    pl = group_index %>% table %>% sqrt %>% as.vector()
    
    #### Initial value
    
    tmp_z = rep(1, 2*p)/p
    tmp_u = rep(0, 2*p+1)
    tmp_beta = rep(1, p)/p
    tmp_beta_0 = 1e-5
    tmp_beta_tilde = c(tmp_beta_0, tmp_beta) 
  }
  
  j = 1
  
  while(j <= iter){
    
    #### Step 1
    i = 1
    loss = Inf
    while(i <= 100){
      
      tmp_W = if_else(abs(((X_tilde %*% tmp_beta_tilde)*4)) <= 1e-8, 1e+8, 
                      abs(1/((X_tilde %*% tmp_beta_tilde)*4))) %>% as.vector() %>% diag()
      
      tmp_beta_tilde = (-1/2) * solve(t(X_tilde) %*% tmp_W %*% X_tilde + (rho/2) * t(Amat %*% K) %*% (Amat %*% K)) %*%
        ((tau - 1/2) * t(X_tilde) %*% rep(1, n) + rho * t(Amat %*% K) %*% ((1/rho) * tmp_u + Bmat %*% tmp_z - Cmat))
      
      loss_new = tau * t(X_tilde %*% tmp_beta_tilde) %*% I(X_tilde %*% tmp_beta_tilde >= 0) - 
        (1-tau) * t(X_tilde %*% tmp_beta_tilde) %*% I(X_tilde %*% tmp_beta_tilde < 0) +
        t(tmp_u) %*% (Amat %*% K %*% tmp_beta_tilde + Bmat %*% tmp_z - Cmat) +
        (rho/2) * norm_vec(Amat %*% K %*% tmp_beta_tilde + Bmat %*% tmp_z - Cmat)^2
      
      if (abs(loss_new - loss) <1e-4) break
      
      loss = loss_new
      
      i = i +1
    }
    
    #### Step 2
    
    tmp_beta = tmp_beta_tilde[-1]
    v = (Amat %*% tmp_beta - Cmat)
    v1 = v[1:p]
    tmp_z1 = v1 + tmp_u[1:p]/rho
    tmp_z1 = if_else(abs(tmp_z1) <= lambda_1/rho, 0, tmp_z1 - lambda_1 * sign(tmp_z1 - lambda_1)/rho)
    
    #### Step 3
    
    v2 = v[(p+1):(2*p)]
    u2 = tmp_u[(p+1):(2*p)]
    
    tmp_z2 = v2 + u2/rho
    
    for (l in unique(group_index)){
      pll = pl[l]
      v2l = v2[group_index == l]
      u2l = u2[group_index == l]
      
      if(norm_vec(tmp_z2[group_index == l]) > lambda_2 * pll/rho){
        tmp_z2[group_index == l] = tmp_z2[group_index == l] -
          lambda_2 * pll * tmp_z2[group_index == l] /(norm_vec(tmp_z2[group_index == l])*rho)
      } else {tmp_z2[group_index == l] = 0 }
    }
    
    tmp_z = c(tmp_z1, tmp_z2)
    
    #### Step 4
    tmp_u = tmp_u + rho * (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)
    
    kkt4 = all(abs(Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) <= 1e-4)
    
    j =  j + 1
    
  }
    
  return(list(current_solution = list(beta = c(tmp_beta_tilde), z = c(tmp_z), u = c(tmp_u))))
}
  
