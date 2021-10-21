library(tidyverse)
library(stats)
library(truncnorm)
library(Matrix)
library(mvtnorm)
library(mvnfast)
library(matrixcalc)
library(LaplacesDemon)
library(purrr)
Rcpp::sourceCpp("real_cpp_fns.cpp")
source("new_mod_fns.R")

## J: number of species in community 
## ng: number of BLOBs/grid cells in region
## TT: number of time points (ex. years) 
## Y: TT-length list of n_{t} x J observation matrices of CPA, where n_{t} are the number of BLOBs with data in t-th time point
## H_ls: (TT-1)-length list of redistribution matrices H_{t}, where the t-th matrix H_{t} is (J * n_{t+1}) x (J * n_{t})
## obs_ids_ls: TT-length list, where t-th entry is a vector of the BLOB ids with data in time point t
## X: TT-length list of n_{t} x p_{x} design matrices for ESI
## Q: TT-length list of n_{t} x p_{q} design matrices for zeros due to unsuitability
## Z: TT-length list of n_{t} x p_{z} design matrices for zeros due to chance
## px, pq, pz: number of predictors in X, Q, Z
## nzA_vec: vector of locations of non-zero elements in A
## nzP_vec: vector of locations of non-zero elements in P

### PRIORS
## w_star_prop_sd_ls: list of lenght TT, where each element contains priors for w_star
## q_eta (r_eta): J-vector of inverse-gamma prior shape (rate) for sigma2_{eta}
## q_gamma (r_gamma): J-vector inverse-gamma prior shape (rate) for sigma2_{gamma}
## A_lb (P_lb): matrix of lower bounds for uniform priors on the alpha (rho) matrix, in format detailed in Supplementary Information
## A_ub (P_ub):  matrix of upper bounds for uniform prior on the alpha (rho) matrix, in format detailed in Supplementary Information


## burnin: number of iterations for burnin
## G: number of iterations post-burnin
## thin_amt: number of iterations to thin by

R <- J*px # number of rows of P matrix; 
X <- lapply(X, t)
z_ids_ls <- lapply(Y, function(x){which(x == 0)})
n_pos <-  apply(do.call(rbind, Y), 2, function(x){sum(x > 0)})


rep_id <- c(); count <- 1
for(j in 1:J){
 rep_id <- c(rep_id, J*count + 1:count)
 count <- count + 1
}

zeroA_vec <- (1:length(A_lb))[-nzA_vec]
A_lb_vec <- A_lb[nzA_vec]; A_ub_vec <- A_ub[nzA_vec]
zeroP_vec <- (1:length(P_lb))[-nzP_vec]
P_lb_vec <- P_lb[nzP_vec]; P_ub_vec <- P_ub[nzP_vec]
U_ids <- (1:J^2)[-rep_id]

update_s2_eta <- T
update_s2_gamma <- T
update_A <- T
update_P <- T
update_delta <- T
update_beta <- T
update_w <- T
update_w_star <- T

## INITIALIZE PARAMETERS THOUGHTFULLY
if(update_s2_eta){
 s2_eta <- r_eta/q_eta 
 s_eta_vec_ls <- list()
 for(t in 1:TT){   
  ng_t <- length(obs_ids_ls[[t]])
  s_eta_vec_ls[[t]] <- sqrt(vec(sapply(s2_eta, function(x){rep(x,ng_t)})))
 }
 print("update s2_eta")
}
if(update_s2_gamma){
 s <- r_gamma/q_gamma 
 s_vec_ls <- list()
 for(t in 1:TT){   
  ng_t <- length(obs_ids_ls[[t]])
  s_vec_ls[[t]] <- sqrt(vec(sapply(s, function(x){rep(x,ng_t)})))
 }
 S <- diag(s)
 print("update s2_gamma")
}
if(update_A){
 A <- matrix(0, nrow = L, ncol = J); 
 A[nzA_vec] <- -0.02
 A[aii_inds] <- 0
 Astart <- A
 Avec <- matrixcalc::vec(A)
 print("update A")
}
if(update_P){
 P <- matrix(0, nrow = R, ncol = J); 
 Pstart <- P
 Pvec <- matrixcalc::vec(P)
 print("update P")
}

if(update_delta){
 delta <- matrix(0, nrow = pq, ncol = J)
 print("update delta")
}
if(update_beta){
 beta <- matrix(0, nrow = pz, ncol = J)
 print("update beta")
}

if(update_w | update_w_star){
 if(update_w_star){
  W_star <- list()
  wstar_ub_ls <- wstar_lb_ls <- list() 
  for(t in 1:(TT-1)){
   ng_t <- length(obs_ids_ls[[t]])
   ub_mat <- matrix( log(apply(Y[[t+1]], 2, max) + 2*s_eta_vec_ls[[t]]),   ncol = J, nrow= ng_t,byrow = T)
   lb_mat_temp <- log(Y[[t+1]])
   lb_mat_temp[which(!is.finite(lb_mat_temp))] <- 100
   lb_mat <-  matrix( apply(lb_mat_temp, 2, min) - 2*s_eta_vec_ls[[t]],   ncol = J, nrow= ng_t,byrow = T) 
   wstar_ub_ls[[t]] <- round(ub_mat,2)
   wstar_lb_ls[[t]] <- round(lb_mat,2)
  }
  print("update w_star")
 }
 if(update_w){
  print("update w_tilde")
  w_ub <- ceiling(log(max(unlist(Y)))) + 2
  w_lb <- 0
  W <- V <- U <- list()
  W_nondeg <- list()
  deg_id_ls <- list()
  for(t in 1:(TT)){
   # randomly choose which are degenerate
   ng_t <- length(obs_ids_ls[[t]])
   draw_deg_mat <- matrix(0,nrow = TT, ncol = ng_t*J)
   deg_temp <- matrix(0, nrow = ng_t, ncol = J)
   z_ids <- z_ids_ls[[t]]
   QA <- pnorm(Q[[t]]%*%delta)
   ZB<- pnorm(Z[[t]] %*% beta)
   prob_deg_cond = (QA)/(QA + (1-QA)*ZB)
   for(i in z_ids){
    if(rbernoulli(1, prob_deg_cond[i])){
     deg_temp[i] <- 1
     draw_deg_mat[t,i] <- 1
    } else{
     deg_temp[i] <- 2
    }
   }
   deg_id_ls[[t]] <- deg_temp
  }
 } 
 
 W_star_deg <- list()
 W_star_nondeg <- list()
 for(t in 1:(TT-1)){
  ng_t <- length(obs_ids_ls[[t]])
  if(update_w){
   if(t== 1){
    W[[t]] <- W_nondeg[[t]] <- matrix(0, nrow = ng, ncol = J)
    W1_init <-  log(Y1_init)
    n_neg <- length(which(W1_init < 0)); n_inf <- length(which(is.infinite(W1_init)))
    W1_init[W1_init<0] <- runif(n_neg, 0, 0.1)
    W1_init[is.infinite(W1_init)] <- runif(n_inf, 0, 0.1)
    W[[t]] <-W_nondeg[[t]] <-W1_init
   } else{
    ng_tm1 <- length(obs_ids_ls[[t-1]])
    ng_t <- length(obs_ids_ls[[t]])
    Wstar_cens <- W_star[[t-1]]
    Wstar_cens[which(Wstar_cens <0)] <- 0
    W[[t]] <- W_nondeg[[t]] <- matrix(H_ls[[t-1]] %*% vec(Wstar_cens), nrow = ng_t, ncol = J)
    
   }
   W[[t]][which(deg_id_ls[[t]] == 1)] <- 0
   U[[t]] <-  apply(W[[t]], 1, function(x){(x %x% x)[-rep_id]})
   V[[t]] <- apply(matrix(1:ng_t, ncol = 1), 1, function(x){W[[t]][x,] %x% X[[t]][,x]})
  }
  if(update_w_star){
   mm <- as.matrix((W[[t]] +  t(V[[t]]) %*% P + t(U[[t]]) %*% A) + rnorm(ng_t*J, 0, sd =  s_vec_ls[[t]]))
   mm[which(mm < wstar_lb_ls[[t]])] <- wstar_lb_ls[[t]][which(mm < wstar_lb_ls[[t]])]; mm[which(mm > wstar_ub_ls[[t]])]<-wstar_ub_ls[[t]][which(mm > wstar_ub_ls[[t]])]
   deg_ids <- deg_id_ls[[t]]
   temp <- matrix(0, nrow = ng_t, ncol = J)
   for(i in 1:(ng_t*J)){
    if(deg_ids[i] == 1){
     temp[i] <- 0
    } else{
     temp[i] <- mm[i]
    }
   }
   W_star[[t]] <- temp
   W_star_nondeg[[t]] <- mm
   W_star_deg[[t]] <- matrix(0, nrow = ng_t, ncol = J)
  }
 }
 if(update_w){
  t <- TT
  ng_t <- length(obs_ids_ls[[t]])
  Wstar_cens <- W_star[[t-1]]
  Wstar_cens[which(Wstar_cens <0)] <- 0
  W[[t]] <- W_nondeg[[t]] <- matrix(H_ls[[t-1]] %*% vec(Wstar_cens), nrow = ng_t, ncol = J)
  
  W[[t]][which(deg_id_ls[[t]] == 1)] <- 0
 }
}


## store results
store <-  0
S2_GAMMA <- matrix(NA,nrow = G/thin_amt,ncol = J )
S2_ETA <-matrix(NA, nrow = G/thin_amt, ncol = J)
ALPHA <- matrix(NA, nrow = G/thin_amt, ncol = pq*J)
BETA <- matrix(NA, nrow = G/thin_amt, ncol = pz*J)
A_POST <- matrix(NA, nrow = G/thin_amt, ncol = length(A))
P_POST <- matrix(NA, nrow = G/thin_amt, ncol = length(P))
W_POST_NONDEG <-  W_POST_NONDEG2 <- W_STAR_NONDEG <- W_STAR_DEG <-W_STAR_NONDEG2 <- W_STAR_DEG2 <-list()
w_mod_amount <- 20
mod_amount <- 200
prop_deg_ls <- list()
for(t in 1:TT){
 ng_t <-  length(obs_ids_ls[[t]])
 W_POST_NONDEG[[t]] <-  W_POST_NONDEG2[[t]] <-  prop_deg_ls[[t]] <- matrix(0, nrow = ng_t, ncol = J)
 if(t < TT){
  W_STAR_NONDEG[[t]] <- W_STAR_DEG[[t]] <-W_STAR_NONDEG2[[t]] <- W_STAR_DEG2[[t]] <- matrix(0, nrow = ng_t, ncol = J)
 }
}



w_star_prop_sd_ls <-list()
acc_nondeg <- nondeg_count <- list() 
for(t in 1:(TT-1)){
 ng_t <-  length(obs_ids_ls[[t]])
 acc_nondeg[[t]] <- nondeg_count[[t]] <-  rep(0, ng_t * J)
 temp <- matrix(0.3, ncol = J, nrow = ng_t)
 temp[z_ids_ls[[t]]] <- 0.5
 w_star_prop_sd_ls[[t]] <- temp
}

### SAMPLER
for(g in 1:(G+burnin)){
 if(g %% mod_amount == 0 & g < G){
  for(t in 1:(TT-1)){
   change_ids <- which(acc_nondeg[[t]]/(w_mod_amount*nondeg_count[[t]]) > 0.2 & nondeg_count[[t]] > 0.1*mod_amount/w_mod_amount)
   w_star_prop_sd_ls[[t]][change_ids] <- w_star_prop_sd_ls[[t]][change_ids]*1.15
   change_ids <- which(acc_nondeg[[t]]/(w_mod_amount*nondeg_count[[t]]) < 0.1 & nondeg_count[[t]] >  0.1*mod_amount/w_mod_amount) 
   w_star_prop_sd_ls[[t]][change_ids] <- w_star_prop_sd_ls[[t]][change_ids]*0.9
  }
  acc_nondeg <- nondeg_count <- list() 
  for(t in 1:(TT-1)){
   ng_t <-  length(obs_ids_ls[[t]])
   acc_nondeg[[t]] <- nondeg_count[[t]] <-  rep(0, ng_t * J)
  }
 }
 if(update_A){
  A <- .update_A(TT, W, W_star, U, V, S, A, P, zeroA_vec, nzA_vec, A_lb_vec, A_ub_vec, L, J)
 }
 
 if(update_P){ 
  P <- .update_P(TT, W, W_star, U, V, S, P, A, zeroP_vec, nzP_vec, P_lb_vec, P_ub_vec, R, J)
 }
 
 if(update_w){
  for(t in 1:(TT-2)){
   ng_t <-  length(obs_ids_ls[[t]])
   ng_tp1 <-  length(obs_ids_ls[[t+1]])
   if(g %% w_mod_amount == 0){
    deg_id_ls[[t]]  <- .update_deg(t, ng_t, J, deg_id_ls, z_ids_ls, Q, Z, delta, beta)
    deg_id_t <-  deg_id_ls[[t]]
    # number of times non deg
    nondeg_count[[t]][which(deg_id_t != 1)] <- nondeg_count[[t]][which(deg_id_t != 1)] + 1
   }
   deg_id_t <-  deg_id_ls[[t]] 
   if( t == 1){
    Wt <- W[[t]]
    redist <- W1_init
    Wt[which(deg_id_t == 1)] <- 0
    Wt[which(deg_id_t == 2)] <- redist[which(deg_id_t == 2)]
    W[[t]] <- Wt
    U[[t]] <-  apply(W[[t]], 1, function(x){(x %x% x)[-rep_id]})
    V[[t]] <- apply(matrix(1:ng_t, ncol = 1), 1, function(x){W[[t]][x,] %x% X[[t]][,x]})
   } else{
    W[[t]] <- matrix(0, nrow = ng_t, ncol = J)
    W[[t]][which(deg_id_t != 1)] <- W_nondeg[[t]][which(deg_id_t != 1)]
    U[[t]] <-  apply(W[[t]], 1, function(x){(x %x% x)[-rep_id]})
    V[[t]] <- apply(matrix(1:ng_t, ncol = 1), 1, function(x){W[[t]][x,] %x% X[[t]][,x]})
   }
   ret <-  sampW_star_miss(W_tilde = W, W_star = W_star, W_star_deg = W_star_deg,
                           W_star_nondeg = W_star_nondeg, deg_id_ls = deg_id_ls, z_ids_ls = z_ids_ls,
                           Y_ls = Y, ng_t = ng_t, ng_tp1 = ng_tp1, J = J, w_prop_sd = w_star_prop_sd_ls[[t]], s_vec_ls = s_vec_ls, 
                           s_eta_vec_ls = s_eta_vec_ls,
                           t = t-1, TT = TT - 1,U_ids = U_ids-1, V_ids = (1:R)-1, 
                           V_ls = lapply(V, t), U_ls = lapply(U,t), XV = X,
                           wstar_lb = wstar_lb_ls[[t]], wstar_ub = wstar_ub_ls[[t]], H_ls= H_ls, 
                           P = P, A = A, order = (1:(ng_t*J))-1, obs_ids_ls = obs_ids_ls, oneH = F)
   W_star[[t]] <- ret[[1]];
   W_star_deg[[t]][which(deg_id_t == 1)] <- ret[[1]][which(deg_id_t == 1)]
   W_star_nondeg[[t]][which(deg_id_t != 1)] <- ret[[1]][which(deg_id_t != 1)]
   # acc of non deg
   acc_nondeg[[t]][which(deg_id_t != 1)] <- acc_nondeg[[t]][which(deg_id_t != 1)] + ret[[5]][which(deg_id_t != 1)] 
   # number of times non deg

   deg_id_tp1 <- deg_id_ls[[t+1]]
   W_nondeg[[t+1]][which(deg_id_tp1 != 1)] <- ret[[2]][which(deg_id_tp1 != 1)]
   W[[t+1]] <- ret[[2]];
   U[[t+1]] <- t(ret[[3]])
   V[[t+1]] <- t(ret[[4]])
  }
  # update [TT-1]
  t <- TT -1
  ng_t <-  length(obs_ids_ls[[t]])
  ng_tp1 <-  length(obs_ids_ls[[t+1]])
  if(g %% w_mod_amount == 0){
   deg_id_ls[[t]]  <- .update_deg(t, ng_t, J, deg_id_ls, z_ids_ls, X, Z, delta, beta)
   deg_id_t <-  deg_id_ls[[t]]
   # number of times non deg
   nondeg_count[[t]][which(deg_id_t != 1)] <- nondeg_count[[t]][which(deg_id_t != 1)] + 1
  }
  deg_id_t <-  deg_id_ls[[t]]
  W[[t]] <- matrix(0, nrow = ng_t, ncol = J)
  W[[t]][which(deg_id_t != 1)] <- W_nondeg[[t]][which(deg_id_t != 1)]
  U[[t]] <-  apply(W[[t]], 1, function(x){(x %x% x)[-rep_id]})
  V[[t]] <- apply(matrix(1:ng_t, ncol = 1), 1, function(x){W[[t]][x,] %x% X[[t]][,x]})
  ret <-   sampW_TTm1_star_miss(W_tilde = W, W_star = W_star, W_star_deg = W_star_deg,
                                W_star_nondeg = W_star_nondeg, deg_id_ls = deg_id_ls, z_ids_ls = z_ids_ls,
                                Y_ls = Y, ng_t = ng_t, ng_tp1 = ng_tp1, J = J, w_prop_sd = w_star_prop_sd_ls[[t]], 
                                s_vec_ls = s_vec_ls, s_eta_vec_ls = s_eta_vec_ls,
                                t = (TT-1)-1, TT = TT - 1,U_ids = U_ids-1, V_ids = (1:R)-1,
                                V_ls = lapply(V, t), U_ls = lapply(U,t), XV = X,
                                wstar_lb = wstar_lb_ls[[t]], wstar_ub = wstar_ub_ls[[t]], H_ls = H_ls,  
                                P = P, A = A,  order = (1:(ng_t*J))-1, oneH = F)
  W_star[[t]] <- ret[[1]];
  W_star_deg[[t]][which(deg_id_t == 1)] <- ret[[1]][which(deg_id_t == 1)]
  W_star_nondeg[[t]][which(deg_id_t != 1)] <- ret[[1]][which(deg_id_t != 1)]
  # acc of non deg
  acc_nondeg[[t]][which(deg_id_t != 1)] <- acc_nondeg[[t]][which(deg_id_t != 1)] + ret[[3]][which(deg_id_t != 1)] 
  
  W[[t+1]] <- ret[[2]];
  deg_id_tp1 <- deg_id_ls[[t+1]]
  W_nondeg[[t+1]][which(deg_id_tp1 != 1)] <- ret[[2]][which(deg_id_tp1 != 1)]
  
  t <- TT
  ng_t <-  length(obs_ids_ls[[t]])
  if(g %% w_mod_amount == 0){
   deg_id_ls[[t]]  <- .update_deg(t, ng_t, J, deg_id_ls, z_ids_ls, Q, Z, delta, beta)
  }
  deg_id_t <-  deg_id_ls[[t]] 
  W[[t]] <- matrix(0, nrow = ng_t, ncol = J)
  W[[t]][which(deg_id_t != 1)] <- W_nondeg[[t]][which(deg_id_t != 1)]
 }
 
 
 if(update_s2_gamma){
  for(j in 1:J){
   qa <- q_gamma[j] + length(unlist(obs_ids_ls))/2
   to_add <- 0
   for(t in 1:(TT-1)){
    to_add <- to_add + crossprod(W_star[[t]][,j] - (W[[t]] + t(V[[t]]) %*% P + t(U[[t]]) %*% A)[,j])
   }
   ra <-r_gamma[j] + to_add/2
   s[j] <-  1/rgamma(1,qa, ra)
  }
  for(t in 1:TT){
   s_vec_ls[[t]] <-  sqrt(vec(sapply(s, function(x){rep(x, length(obs_ids_ls[[t]]))})))
  }
  S <- diag(s)
 }
 
 if(update_s2_eta){
  # update s2_eta
  for(j in 1:J){
   to_add <- 0
   for(t in 1:TT){
    to_add <- to_add + crossprod((log(Y[[t]][,j])- W[[t]][,j])[which(deg_id_ls[[t]][,j] ==0)])
   }
   s2_eta[j] <- 1/rgamma(1, q_eta[j] + 0.5*n_pos[j], r_eta[j] + 0.5*to_add)
   for(t in 1:TT){
    s_eta_vec_ls[[t]] <-  sqrt(vec(sapply(s2_eta, function(x){rep(x, length(obs_ids_ls[[t]]))})))
   }
  }
 }
 
 if(update_delta){
  ## update delta params
  delta <- .update_delta(pq, J, TT, Q, Z, delta, beta, Y, W, s_eta_vec_ls, deg_id_ls,
                         a0_prior = matrix(rep(-1,J), nrow = 1))
 }
 
 if(update_beta){
  beta <- .update_beta(pz, J, TT, Q, Z, delta, beta, Y, W, s_eta_vec_ls, deg_id_ls)
 }
 
 ## STORE
 if(g > burnin){
  if (g%%thin_amt == 0){
   store <- store + 1
   A_POST[store,] <- A; P_POST[store,] <- P
   ALPHA[store,] <- delta; BETA[store,] <- beta
   S2_ETA[store,] <- s2_eta; S2_GAMMA[store,] <- s
  }
   for(t in 1:TT){
    W_POST_NONDEG[[t]] <- W_POST_NONDEG[[t]] + W_nondeg[[t]]
    W_POST_NONDEG2[[t]] <- W_POST_NONDEG2[[t]] + (W_nondeg[[t]])^2
    prop_deg_ls[[t]][which(deg_id_ls[[t]] == 1)] <- prop_deg_ls[[t]][which(deg_id_ls[[t]] == 1)] + 1
    if(t < TT){
     W_STAR_NONDEG[[t]] <- W_STAR_NONDEG[[t]] + W_star_nondeg[[t]]
     W_STAR_NONDEG2[[t]] <- W_STAR_NONDEG2[[t]] + (W_star_nondeg[[t]])^2
     W_STAR_DEG[[t]] <- W_STAR_DEG[[t]] + W_star_deg[[t]]
     W_STAR_DEG2[[t]] <- W_STAR_DEG2[[t]] + (W_star_deg[[t]])^2
    }
   }
 
 }
} 

