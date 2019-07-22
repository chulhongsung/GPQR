rm(list = ls())
gc(reset = T)

library(dplyr)
setwd("~/Github/SGL")
load("./data/example_kospi_dat.RData")
norm_vec <- function(x) sqrt(sum(x^2))


p = 20
n = 100
# convert the raw data into log return rate
return_table = lapply(example_kospi_dat, 
                      function(x) log(lag(x, 1)/lag(x, 0))) %>% 
               bind_cols() %>% as.matrix() %>% .[-1,1:p] %>% .[1:n,]
data = return_table*100

# model parameter in the objective function
tau = 0.5
lambda_1 = 5
lambda_2 = 0

# constraint (mu_zero:target return rate, group_index: group index)
mu = colMeans(data)
mu_zero = max(mu)
group_index = c(rep(1,8), rep(2,12))
pl = group_index %>% table %>% sqrt %>% as.vector()

# ADMM optimization parameters
rho = n
iter_outer = 1e+6
inner_iter = 100


# ADMM input
Amat = rbind(diag(p), diag(p), 1, c(mu))
Bmat = rbind(-diag(2*p), 0, 0)
Cmat = c(rep(0, 2*p), 1, mu_zero)
X_tilde = cbind(-1, data)
K = cbind(0, diag(p)) 

#### Initializing parameter in ADMM
tmp_z = rep(0, 2*p)/p
tmp_u = rep(0, 2*p+2)
tmp_beta = rep(1, p)/p
tmp_beta_0 = 1e-5
tmp_beta_tilde = c(tmp_beta_0, tmp_beta)

j = 1
while(j <= iter_outer)
{
  #### Loss function
#  tmp_z1 = tmp_z[1:p]
#  tmp_z2 = tmp_z[(p+1):(2*p)]
  #### Step 1
  i = 1
  loss = Inf
  while(i <= inner_iter)
  {
    tmp_W = if_else(abs(((X_tilde %*% tmp_beta_tilde)*4)) <= 1e-8, 1e+8, 
                    abs(1/((X_tilde %*% tmp_beta_tilde)*4))) %>% as.vector() %>% diag()
    
    tmp_beta_tilde = (-1/2) * solve(t(X_tilde) %*% tmp_W %*% X_tilde + 
                                      (rho/2) * t(Amat %*% K) %*% (Amat %*% K)) %*%
                              ((tau - 1/2) * t(X_tilde) %*% rep(1, n) + 
                                 rho * t(Amat %*% K) %*% ((1/rho) * tmp_u + Bmat %*% tmp_z - Cmat))
   # loss_for_beta0 = function(X, tmp_beta, tmp_beta_zero, tau, tmp_u, tmp_z, Amat, Bmat, Cmat){
   #    tau * t(X %*% tmp_beta - tmp_beta_zero) %*% I(X %*% tmp_beta - tmp_beta_zero >= 0) - (1-tau) * t(X %*% tmp_beta - tmp_beta_zero) %*% I(X %*% tmp_beta - tmp_beta_zero < 0) +
   #      t(tmp_u) %*% (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) +
   #      (rho/2) * norm_vec(Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)^2 }
   # 
   #  tmp_beta_zero = optimize(loss_for_beta0, c(-10, 10), X = return_table, tmp_beta = tmp_beta_tilde[-1], tau = tau, tmp_z = tmp_z, tmp_u = tmp_u, Amat = Amat, Bmat = Bmat, Cmat = Cmat)$minimum
   # 
   #  tmp_beta_tilde = c(tmp_beta_zero, tmp_beta_tilde[-1])
   #  
    loss_new = tau * t(X_tilde %*% tmp_beta_tilde) %*% I(X_tilde %*% tmp_beta_tilde >= 0) - (1-tau) * t(X_tilde %*% tmp_beta_tilde) %*% I(X_tilde %*% tmp_beta_tilde < 0) +
      t(tmp_u) %*% (Amat %*% K %*% tmp_beta_tilde + Bmat %*% tmp_z - Cmat) +
      (rho/2) * norm_vec(Amat %*% K %*% tmp_beta_tilde + Bmat %*% tmp_z - Cmat)^2
    if (abs(loss_new - loss) <1e-4) break
    loss = loss_new
    #print(loss_)
    # gradient_beta = tau * t(X_tilde) %*% I(X_tilde %*% tmp_beta_tilde >= 0) - (1 - tau) * t(X_tilde) %*% I(X_tilde %*% tmp_beta_tilde < 0) +
    #   rho * t(Amat %*% K) %*% (Amat %*% K) %*% tmp_beta_tilde + rho * t(Amat %*% K) %*% ((1/rho) * tmp_u + Bmat %*% tmp_z - Cmat)
    #print(gradient_beta/n)
    i = i +1
  }
  
  #### Step 2
  # gradient_beta = tau * t(X_tilde) %*% I(X_tilde %*% tmp_beta_tilde >= 0) - (1 - tau) * t(X_tilde) %*% I(X_tilde %*% tmp_beta_tilde < 0) + 
  #   rho * t(Amat %*% K) %*% (Amat %*% K) %*% tmp_beta_tilde + rho * t(Amat %*% K) %*% ((1/rho) * tmp_u + Bmat %*% tmp_z - Cmat)
  
  #kkt1 = all(abs(gradient_beta) <= 1e-4)
  tmp_beta = tmp_beta_tilde[-1]
  v = (Amat %*% tmp_beta - Cmat)
  v1 = v[1:p]
  tmp_z1 = v1 + tmp_u[1:p]/rho
  tmp_z1 = if_else(abs(tmp_z1) <= lambda_1, 0, tmp_z1 - lambda_1 * sign(tmp_z1 - lambda_1))
  #kkt2 = all(abs(rho * (tmp_z1 - v1) - tmp_u[1:p]) <= lambda_1)
  
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
  
  #kkt3 = all(tapply(tmp_z2, group_index, function(x) norm_vec(x)) <= lambda_2)
  
  tmp_z = c(tmp_z1, tmp_z2)
  
  #### Step 4
  tmp_u = tmp_u + rho * (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)
  aa = Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat
  #kkt4 = all(abs(Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) <= 1e-4)
  if (j%%(1e+3) == 0) 
  {
    cat(j,'\n')
    print(tail(aa,2))
    cat("coef::",round(tmp_beta_tilde[-1],3) ,'\n')
    loss = sum(tau * (X_tilde %*% tmp_beta_tilde) *I(X_tilde %*% tmp_beta_tilde >= 0) + 
                   (1-tau)*(-X_tilde %*% tmp_beta_tilde)*I(X_tilde %*% tmp_beta_tilde < 0)) + 
        lambda_1 * sum(abs(tmp_z1)) + 
        lambda_2 * sum(pl * tapply(tmp_z2, group_index, function(x) norm_vec(x))) + 
        t(tmp_u) %*% (Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat) + 
        (rho/2) * norm_vec(Amat %*% tmp_beta + Bmat %*% tmp_z - Cmat)^2
    print(loss)
  }
  # if (j%%(1e+4) == 0)
  # {
  #   rho = 2*rho
  # }
  j =  j + 1
}

sum(tmp_beta_tilde[-1])
tmp_beta_tilde
sum(tmp_beta_tilde[-1]*mu)
mu_zero

