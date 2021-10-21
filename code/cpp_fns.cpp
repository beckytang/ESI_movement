#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
int mod(int x, int y){
 return x - floor(x/y)*y;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat vec_cpp(arma::mat W0){
  return W0.as_col();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat vec2mat(arma::vec x, int nrow, int ncol) {
  arma::mat y(x);
  y.reshape(nrow, ncol);
  return y;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
bool in_cpp(arma::mat x, double y){
 bool ret = true;
 arma::uvec idx = find(x == y);
 double len = idx.size();
 if(len == 0){
  ret = false;
 }
 return ret;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec tapply_sum_cpp(arma::vec x, arma::vec x_index, int ng2){
  arma::vec ret = zeros(ng2);
  arma::uvec idx;
  for(int i = 0; i < ng2; i++){
   idx = find(x_index == i);
    ret[i] = sum(x.elem(idx));
  }
  return(ret);
}

// [[Rcpp::depends(RcppArmadillo)]]
double tnormRcpp(double lo, double hi, double mu, double sig){
 
 double q1, q2, z;
 
 q1 = Rf_pnorm5(lo,mu,sig,1,0);
 q2 = Rf_pnorm5(hi,mu,sig,1,0);
 z = q1 + unif_rand()*(q2-q1);
 z = Rf_qnorm5(z, mu, sig, 1, 0);
 
 if(z > hi){
  z = lo;
 }
 
 if(z < lo){
  z = hi;
 }
 return(z);
}

// [[Rcpp::export]]
arma::rowvec kron(arma::vec A, arma::vec B, arma::uvec ret_ids){
 int n = A.n_elem;
 int m = B.n_elem;
 int r = ret_ids.n_elem;
 arma::vec ret(n*m);
 arma::rowvec row_ret(r);
 int count=0;
 for(int i = 0; i < n; i++){
  for(int j = 0; j < m; j++){
   ret[count] = A[i]*B[j];
   count++;
  }
 }
 row_ret = ret.elem(ret_ids).t();
 return row_ret;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat rmvnormRcpp(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  bool success = false;
  arma::mat S = sigma;
  int i = 1;
  
  arma::mat Y = randn(n, ncols);
  
  while(success == false && i < 5){
    
    success = chol(S, sigma);
    
    if(success == false){
      sigma += eye(ncols,ncols) * 1e-5;
    }
    
    i = i + 1;
  }
  
  if(success == false){
    //    throw std::range_error("sigma not positive definite");
    return arma::repmat(mu*0, 1, n).t();
  }
  
  return arma::repmat(mu, 1, n).t() + Y * chol(sigma);
}



// [[Rcpp::export]]
double dtnorm(double x, double lo, double hi, double mu, double sigma){
 double q1, q2, ret=0;
 q1 = Rf_pnorm5(lo,mu,sigma,1,0);
 q2 = Rf_pnorm5(hi,mu,sigma,1,0);
 if (x > lo & x < hi){
  ret =  R::dnorm(x, mu, sigma, FALSE)/(q2-q1);
 }
 return ret;
 
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_star(const List& W_tilde, const List& W_star, const List& W_star_deg,
                const List& W_star_nondeg, const List& deg_id_ls, 
                 const List& z_ids_ls, const List& Y_ls, arma::mat W_star_true_t,
                 int ng, int J, arma::mat w_prop_sd, arma::vec s_vec,
                 double s2_eta, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                 arma::mat W0_tilde, 
                 const List& V_ls, const List& U_ls, const List& XV,
                 arma::mat wstar_ub, arma::mat wstar_lb,
                 sp_mat H, arma::mat P, arma::mat A, 
                 IntegerVector order, bool oneH = true) {
  IntegerVector ord = order;
  arma::vec all_props = zeros(ng*J);
  arma::vec LH1 = zeros(ng*J), LH2 = zeros(ng*J), LH3= zeros(ng*J), LH4= zeros(ng*J);
  arma::mat Ytp1 = Y_ls[t+1];
  arma::mat W_tp1= W_tilde[t+1];
  arma::mat deg_tp1 = deg_id_ls[t+1];
  arma::uvec deg_id_tp1 = find( deg_tp1 == 1);
  arma::mat deg_t = deg_id_ls[t];
  arma::vec z_tp1_ids = z_ids_ls[t+1];
  arma::mat Utp1 = U_ls[t+1];
  arma::mat Vtp1 = V_ls[t+1];
  arma::mat Wstar_tp1 = W_star[t+1];

  arma::mat XVtp1 = XV[t+1];

  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat Wstar_t = W_star[t];
  arma::mat Wstar_nondeg_t = W_star_nondeg[t];
  arma::mat Wstar_deg_t = W_star_deg[t];
  arma::mat deg_temp = deg_id_ls[t];
  arma::vec acc = zeros(ng*J);
  
  arma::mat Ut = U_ls[t];
  arma::mat Vt = V_ls[t];
  arma::mat XVt = XV[t];
  arma::mat m_t = Wt + Vt* P + Ut*A, mprop_tp1, m_tp1;
  int i, j;
  double lh1 = 0, lh2 = 0, lh3 = 0, lh4 = 0, lhr = 0, u;
  bool is_zero;
  arma::uvec neg_ids1, neg_ids2;
  double wstar_prop_t, wstar_t;
  arma::mat Wstar_prop_t, Wprop_tp1, Wstar_cens_prop, Uprop_tp1 = Utp1, Vprop_tp1 = Vtp1;
  arma::mat Wstar_cens;
  // neg_ids2 = find(Wstar_cens < 0);
  // Wstar_cens.elem(neg_ids2).fill(0);
  // W_tp1 = H * Wstar_cens;
  // W_tp1.elem(deg_id_tp1).fill(0);
  // for(int nn = 0; nn < ng; nn++){
  //   Utp1.row(nn) = kron(W_tp1.row(nn).t(), W_tp1.row(nn).t(), U_ids);
  //   Vtp1.row(nn) = kron(W_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
  // }
  std::random_shuffle(ord.begin(), ord.end());
  for(int q = 0; q < ng*J; q++){
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    lh4 = 0;
    i = ord[q];

    // k = mod(i+1, ng);
    // if(k == 0){
    //   k = ng ;
    // }
    // k = k-1;
    j = ceil((i+1)/ng) - 1;
    
    if(deg_t[i] == 1){
      // when W_tilde = deg. 0
      wstar_t = Wstar_deg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = tnormRcpp(-0.3, 0.3, wstar_t, 0.075);
      lh4 = log(dtnorm(wstar_t, -0.3, 0.3, wstar_prop_t, 0.075))-
        log(dtnorm(wstar_prop_t, -0.3, 0.3, wstar_t, 0.075));
    } else{
      wstar_t = Wstar_nondeg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = R::rnorm(wstar_t, w_prop_sd[i]);
      lh4 = 0;
    }
    if(wstar_prop_t <= wstar_ub[i] & wstar_prop_t >= wstar_lb[i]){
      Wstar_prop_t = Wstar_t;
      Wstar_prop_t[i] = wstar_prop_t;
      // growth
      lh1 = R::dnorm(wstar_prop_t, m_t[i], s_vec[i], TRUE)-
        R::dnorm(wstar_t, m_t[i], s_vec[i], TRUE);
      
      // project forward
      Wstar_cens_prop = Wstar_prop_t;
      neg_ids1 = find(Wstar_cens_prop < 0);
      Wstar_cens_prop.elem(neg_ids1).fill(0);
      if(oneH){
        Wprop_tp1 = H * Wstar_cens_prop;
      } else{
        Wprop_tp1 = vec2mat(H * vec_cpp(Wstar_cens_prop), ng, J);
      }
      Wprop_tp1.elem(deg_id_tp1).fill(0);
      
      Wstar_cens = Wstar_t;
      neg_ids2 = find(Wstar_cens < 0);
      Wstar_cens.elem(neg_ids2).fill(0);
      if(oneH){
        W_tp1 = H * Wstar_cens;
      } else{
        W_tp1 = vec2mat(H * vec_cpp(Wstar_cens), ng, J);
      }
      W_tp1.elem(deg_id_tp1).fill(0);
      // 
      // Y_tp1 ~ LN(W_tile_tp1)
      for(int n = 0; n < ng*J; n++){
        is_zero = in_cpp(z_tp1_ids,n+1); 
        if(is_zero == false){
          lh2 = lh2 + R::dlnorm(Ytp1[n], Wprop_tp1[n], sqrt(s2_eta), TRUE)-
            R::dlnorm(Ytp1[n], W_tp1[n], sqrt(s2_eta), TRUE);
        }
      }
      
      for(int nn = 0; nn < ng; nn++){
        Utp1.row(nn) = kron(W_tp1.row(nn).t(), W_tp1.row(nn).t(), U_ids);
        Vtp1.row(nn) = kron(W_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
        Uprop_tp1.row(nn) = kron(Wprop_tp1.row(nn).t(), Wprop_tp1.row(nn).t(), U_ids);
        Vprop_tp1.row(nn) = kron(Wprop_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
      }
      
      mprop_tp1 =  Wprop_tp1 + Vprop_tp1 * P + Uprop_tp1 * A;
      m_tp1 =  W_tp1 + Vtp1* P + Utp1*A;
      for(int l = 0; l < ng*J; l++){
        lh3 = lh3 + R::dnorm(Wstar_tp1[l], mprop_tp1[l], s_vec[l], TRUE)-
          R::dnorm(Wstar_tp1[l], m_tp1[l], s_vec[l], TRUE);
      }
      
      lhr = lh1 +lh2 + lh3 + lh4;// + 
      // R::dnorm(wstar_prop_t,W_star_true_t[i], 0.05, TRUE)-
      // R::dnorm(wstar_t,W_star_true_t[i], 0.05, TRUE);
      u = runif(1)[0];
      if(lhr > log(u)){
        Wstar_t[i] = wstar_prop_t;
        W_tp1 = Wprop_tp1;
        Utp1 = Uprop_tp1;
        Vtp1 = Vprop_tp1;
        acc[i] = 1;
      }
    }
    LH1(i) = lh1;
    LH2(i) = lh2;
    LH3(i) = lh3;
    LH4(i) = lh4;
    all_props(i) = wstar_prop_t;
  }
  return List::create(Wstar_t, W_tp1, Utp1, Vtp1, acc, LH1, LH2, LH3, LH4, all_props,
                      m_t);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_TTm1_star(const List& W_tilde, const List& W_star, const List& W_star_deg,
                     const List& W_star_nondeg, const List& deg_id_ls, 
                const List& z_ids_ls, const List& Y_ls, arma::mat W_star_true_t,
                int ng, int J, arma::mat w_prop_sd, arma::vec s_vec,
                double s2_eta, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                arma::mat W0_tilde,
                const List& V_ls, const List& U_ls, const List& XV,
                arma::mat wstar_ub, arma::mat wstar_lb,
                sp_mat H, arma::mat P, arma::mat A, 
                IntegerVector order, bool oneH = true) {
  IntegerVector ord = order;
  arma::mat Yt = Y_ls[t];
  arma::mat Ytp1 = Y_ls[t+1];
  arma::mat W_tp1 = W_tilde[t+1];
  arma::mat deg_tp1 = deg_id_ls[t+1];
  arma::uvec deg_id_tp1 = find( deg_tp1 == 1);
  arma::mat deg_t = deg_id_ls[t];
  arma::vec z_tp1_ids = z_ids_ls[t+1];

  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat Wstar_t = W_star[t];
  arma::mat Wstar_deg_t = W_star_deg[t];
  arma::mat Wstar_nondeg_t = W_star_nondeg[t];
  arma::mat deg_temp = deg_id_ls[t];
  arma::vec acc = zeros(ng*J);
  // 
  arma::mat Ut = U_ls[t];
  arma::mat Vt = V_ls[t];
  arma::mat XVt = XV[t];
  arma::mat m_t = Wt + Vt* P + Ut*A;
  int i, j;
  arma::vec LH1 = zeros(ng*J), LH2 = zeros(ng*J), LH3= zeros(ng*J), LH4= zeros(ng*J),
    all_props = zeros(ng*J);
  double lh1 = 0, lh2 = 0, lh3 = 0, lh4 = 0, lhr = 0, u;
  bool is_zero;

  double wstar_prop_t, wstar_t;
  arma::mat Wstar_prop_t, Wprop_tp1, Wstar_cens_prop;

  arma::mat Wstar_cens;
  arma::uvec neg_ids1, neg_ids2 ;
  std::random_shuffle(ord.begin(), ord.end());
  for(int q = 0; q < ng*J; q++){
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    lh4 = 0;
    i = ord[q];
   
    j = ceil((i+1)/ng) - 1;
    
    if(deg_t[i] == 1){
      // Wtilde = deg 0
      wstar_t = Wstar_deg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = tnormRcpp(-0.3, 0.3, wstar_t, 0.075);
      lh4 = log(dtnorm(wstar_t, -0.3, 0.3, wstar_prop_t, 0.075))-
        log(dtnorm(wstar_prop_t, -0.3, 0.3, wstar_t, 0.075));
    } else{
      wstar_t = Wstar_nondeg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = R::rnorm(wstar_t, w_prop_sd[i]);
    }
    if(wstar_prop_t <= wstar_ub[i] & wstar_prop_t >= wstar_lb[i]){
      Wstar_prop_t = Wstar_t;
      Wstar_prop_t[i] = wstar_prop_t;
      
      lh1 = R::dnorm(wstar_prop_t, m_t[i], s_vec[i], TRUE)-
        R::dnorm(wstar_t, m_t[i], s_vec[i], TRUE);
      Wstar_cens_prop = Wstar_prop_t;
      neg_ids1 = find(Wstar_cens_prop < 0);
      Wstar_cens_prop.elem(neg_ids1).fill(0);
      if(oneH){
        Wprop_tp1 = H * Wstar_cens_prop;
      } else{
        Wprop_tp1 = vec2mat(H * vec_cpp(Wstar_cens_prop),ng,J);
      }
      Wprop_tp1.elem(deg_id_tp1).fill(0);
      
      Wstar_cens = Wstar_t;
      neg_ids2 = find(Wstar_cens < 0);
      Wstar_cens.elem(neg_ids2).fill(0);
      if(oneH){
        W_tp1 = H * Wstar_cens;
      } else{
        W_tp1 = vec2mat(H * vec_cpp(Wstar_cens),ng, J);
      }
      W_tp1.elem(deg_id_tp1).fill(0);
      //
      for(int n = 0; n < ng*J; n++){
        is_zero = in_cpp(z_tp1_ids,n+1);
        if(is_zero == false){
          lh2 = lh2 + R::dlnorm(Ytp1[n], Wprop_tp1[n], sqrt(s2_eta), TRUE)-
            R::dlnorm(Ytp1[n], W_tp1[n], sqrt(s2_eta), TRUE);
        }
      }
      
      lhr = lh1 +lh2 + lh4 ;//+
      // R::dnorm(wstar_prop_t,W_star_true_t[i], 0.075, TRUE)-
      // R::dnorm(wstar_t,W_star_true_t[i], 0.075, TRUE);
      u = runif(1)[0];
      if(lhr > log(u)){
        Wstar_t[i] = wstar_prop_t;
        W_tp1 = Wprop_tp1;
        acc[i] = 1;
      } 
    }
    LH1(i) = lh1;
    LH2(i) = lh2;
    all_props(i) = wstar_prop_t;
  }
  return List::create(Wstar_t, W_tp1, acc, LH1, LH2, all_props);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_tilde2(const List& W_tilde, const List& W_star, const List& deg_id_ls,
                 const List& z_ids_ls, const List& Y_ls, const List& W_tilde_nondeg,
                 const List& W_star_deg, const List& W_star_nondeg,
                 int ng, int J, double w_prop_sd, arma::vec s_vec,
                 double s2_eta, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                 arma::mat W0_tilde,
                 const List& V_ls, const List& U_ls, const List& XV, 
                 double w_ub, double w_lb,
                 sp_mat H, arma::mat P, arma::mat A, 
                 IntegerVector order, bool oneH = true) {
  IntegerVector ord = order;
  arma::mat Yt = Y_ls[t];
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat W_star_t = W_star[t];
  arma::mat deg_temp = deg_id_ls[t];
  arma::mat deg_tp1 = deg_id_ls[t+1];
  arma::uvec deg_ids = find(deg_temp == 1);
  arma::vec z_ids = z_ids_ls[t], acc = zeros(ng*J);
  arma::mat Ut = U_ls[t], Uprop_t = Ut;
  arma::mat Vt = V_ls[t], Vprop_t = Vt;
  arma::mat XVt = XV[t];
  arma::mat Wstar_tp1 = W_star[t];
  arma::mat redist = Wt;
  arma::mat Wstar_tm1_cens = W_star[t];
  arma::mat W_nondeg_t = W_tilde_nondeg[t];
  
  int i, j;
  double w_tilde_prop = -1;
  double lh1 = 0, lh2 = 0, lh3 = 0, lhr = 0, u;
  if(t <= (TT-2)){
    arma::mat Wstar_tp1 = W_star[t+1];
  }

  arma::mat W_star_deg_t = W_star_deg[t];
  arma::mat W_star_nondeg_t = W_star_nondeg[t];
  arma::mat mean_curr = Wt + Vt* P + Ut*A;
  arma::mat mean_prop;

  for(int q = 0; q < ng*J; q++){
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    i = ord[q];
    j = ceil((i+1)/ng) - 1;
    
    if(deg_temp[i] == 1){
      w_tilde_prop = 0;
      lh1 = 0;
      lh3 = 0;
    }
    if(deg_temp[i] == 2){
      Wt[i] = W_nondeg_t[i];
      w_tilde_prop = tnormRcpp(0, R_PosInf, Wt[i], w_prop_sd);
      // lh1 = R::dlnorm(Yt[i], w_tilde_prop, sqrt(s2_eta), TRUE) -
      //   R::dlnorm(Yt[i], Wt[i], sqrt(s2_eta), TRUE);
      lh1 = 0;
      lh3 = log(dtnorm(Wt[i], 0 ,  R_PosInf, w_tilde_prop, w_prop_sd))-
        log(dtnorm(w_tilde_prop, 0, R_PosInf, Wt[i], w_prop_sd));
    }
    if(deg_temp[i] == 0){
      Wt[i] = W_nondeg_t[i];
      w_tilde_prop = tnormRcpp(0, R_PosInf, Wt[i], w_prop_sd);
      lh1 = R::dlnorm(Yt[i], w_tilde_prop, sqrt(s2_eta), TRUE) -
        R::dlnorm(Yt[i], Wt[i], sqrt(s2_eta), TRUE);
      lh3 = log(dtnorm(Wt[i], 0 ,  R_PosInf, w_tilde_prop, w_prop_sd))-
        log(dtnorm(w_tilde_prop, 0, R_PosInf, Wt[i], w_prop_sd));
    }
    W_tilde_prop = Wt;
    W_tilde_prop[i] = w_tilde_prop;

    for(int nn = 0; nn < ng; nn++){
      Ut.row(nn) = kron(Wt.row(nn).t(), Wt.row(nn).t(), U_ids);
      Vt.row(nn) = kron(Wt.row(nn).t(), XVt.col(nn), V_ids);
      Uprop_t.row(nn) = kron(W_tilde_prop.row(nn).t(), W_tilde_prop.row(nn).t(), U_ids);
      Vprop_t.row(nn) = kron(W_tilde_prop.row(nn).t(), XVt.col(nn), V_ids);
    }
    mean_curr = Wt + Vt * P + Ut * A;
    mean_prop = W_tilde_prop + Vprop_t * P + Uprop_t*A;
    for(int n = 0; n < ng*J; n++){
      lh2 = lh2+
        R::dnorm( W_star_t[n], mean_prop[n], s_vec[n], TRUE)-
        R::dnorm( W_star_t[n], mean_curr[n], s_vec[n], TRUE);
    } 
    
    if(deg_temp[i] == 1){
      lhr = 1000; // automatic accept
    } else{
      lhr = lh1 + lh2 + lh3;
    }
    
    u = runif(1)[0];
    if(w_tilde_prop >= w_lb & w_tilde_prop <= w_ub){
      if(lhr > log(u)){
        Wt[i] = w_tilde_prop;
        Ut = Uprop_t;
        Vt = Vprop_t;
        acc[i] = 1;
      }
    }

  }
  return List::create(Wt, deg_temp, Ut, Vt, acc);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_t_tilde2(const List& W_tilde, const List& W_star, const List& deg_id_ls,
                   const List& z_ids_ls, const List& Y_ls, const List& W_tilde_nondeg,
                   int ng, int J, double w_prop_sd, arma::vec s_vec,
                   double s2_eta, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                   arma::mat W0_tilde,
                   const List& V_ls, const List& U_ls, const List& XV, 
                   double w_lb, double w_ub,
                   sp_mat H, arma::mat P, arma::mat A,
                   IntegerVector order, bool oneH = true) {
  IntegerVector ord = order;
  arma::mat Yt = Y_ls[t];
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat deg_temp = deg_id_ls[t];
  arma::vec z_ids = z_ids_ls[t];
  arma::mat W_nondeg_t = W_tilde_nondeg[t];
  int i, j;
  double w_tilde_prop;
  double lh1 = 0, lh2 = 0, lhr = 0, u;
  arma::vec acc = zeros(ng*J);
  arma::mat Wstar_tm1_cens = W_star[t-1];
  arma::uvec neg_ids = find(Wstar_tm1_cens < 0);
  Wstar_tm1_cens.elem(neg_ids).fill(0);
  arma::mat redist = H * Wstar_tm1_cens;

  for(int q = 0; q < ng*J; q++){
    i = ord[q];
    j = ceil((i+1)/ng) - 1;

    if(deg_temp[i] == 1){
      w_tilde_prop = 0;
      deg_temp[i] = 1;
      lhr = 1000; //automatic accept
    }
    if(deg_temp[i] == 2){
      Wt[i] = W_nondeg_t[i];
      w_tilde_prop = tnormRcpp(0, R_PosInf, Wt[i], w_prop_sd);
      W_tilde_prop = Wt;
      deg_temp[i] = 2;
      lh1 = 0;
      // lh1 = R::dlnorm(Yt[i], w_tilde_prop, sqrt(s2_eta), TRUE) -
      //   R::dlnorm(Yt[i], Wt[i], sqrt(s2_eta), TRUE);
      lh2 = log(dtnorm(Wt[i], 0 ,  R_PosInf, w_tilde_prop, w_prop_sd))-
        log(dtnorm(w_tilde_prop, 0, R_PosInf, Wt[i], w_prop_sd));    
      lhr = lh1 + lh2;
    }
    if(deg_temp[i] == 0){
      Wt[i] =  W_nondeg_t[i];
      w_tilde_prop = tnormRcpp(0, R_PosInf, Wt[i], w_prop_sd);
      W_tilde_prop = Wt;
      W_tilde_prop[i] = w_tilde_prop;
      lh1 = R::dlnorm(Yt[i], w_tilde_prop, sqrt(s2_eta), TRUE) -
        R::dlnorm(Yt[i], Wt[i], sqrt(s2_eta), TRUE);
      lh2 = log(dtnorm(Wt[i], 0 ,  R_PosInf, w_tilde_prop, w_prop_sd))-
        log(dtnorm(w_tilde_prop, 0, R_PosInf, Wt[i], w_prop_sd));
      lhr = lh1 + lh2;
    }
    u = runif(1)[0];
    if(w_tilde_prop >= w_lb & w_tilde_prop <= w_ub){
      if(lhr > log(u) ){
        Wt[i] = w_tilde_prop;
        acc[i] = 1;
      }
    }
  }

  return List::create(Wt, deg_temp, acc);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_tilde3(const List& W_tilde, const List& W_star, const List& deg_id_ls,
                  const List& z_ids_ls, const List& Y_ls, const List& W_tilde_nondeg,
                  const List& W_star_deg, const List& W_star_nondeg,
                  int ng, int J, double w_prop_sd, arma::vec s_vec,
                  double s2_eta, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                  arma::mat W0_tilde, arma::mat draw_deg_mat, int g,
                  const List& V_ls, const List& U_ls, const List& XV, 
                  double w_ub, double w_lb,
                  sp_mat H, arma::mat P, arma::mat A, 
                  IntegerVector order, bool oneH = true) {
  IntegerVector ord = order;
  arma::mat Yt = Y_ls[t];
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat W_star_t = W_star[t];
  arma::mat deg_temp = deg_id_ls[t];
  arma::mat deg_tp1 = deg_id_ls[t+1];
  arma::uvec deg_ids = find(deg_temp == 1);
  arma::vec z_ids = z_ids_ls[t], acc = zeros(ng*J);
  arma::mat Ut = U_ls[t], Uprop_t = Ut;
  arma::mat Vt = V_ls[t], Vprop_t = Vt;
  arma::mat XVt = XV[t];
  arma::mat Wstar_tp1 = W_star[t];
  arma::mat redist_deg, redist_nondeg, Wstar_tm1_cens_temp;
  NumericVector v(2);

  arma::mat Wstar_tm1_cens = W_star[t-1];
  arma::mat W_star_deg_tm1 = W_star_deg[t-1];
  arma::mat W_star_nondeg_tm1 = W_star_nondeg[t-1];
  arma::uvec neg_ids = find(Wstar_tm1_cens < 0);
  Wstar_tm1_cens.elem(neg_ids).fill(0);
  

  arma::mat W_nondeg_t = W_tilde_nondeg[t];
  //double prop_deg;
  
  int i, j;
  double w_tilde_prop = -1;
  double lh1 = 0, lh2 = 0, lh3 = 0, lhr = 0, u;
  if(t <= (TT-2)){
    arma::mat Wstar_tp1 = W_star[t+1];
  }
  
  arma::mat W_star_deg_t = W_star_deg[t];
  arma::mat W_star_nondeg_t = W_star_nondeg[t];
  
  arma::mat mean_curr, mean_prop;

  for(int q = 0; q < ng*J; q++){
    i = ord[q];
   
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    j = ceil((i+1)/ng) - 1;
    
    if(deg_temp[i] == 1){
      w_tilde_prop = 0;
      lh1 = 0;
      lh3 = 0;
    }
    if(deg_temp[i] == 2){
      // prop_deg = draw_deg_mat[i]/g;
      // if(t > 0){
      //   v = {0, W_star_deg_tm1[i]};
      //   Wstar_tm1_cens_temp = Wstar_tm1_cens;
      //   Wstar_tm1_cens_temp[i] = max(v);
      //   redist_deg = H*Wstar_tm1_cens_temp;
      // 
      //   v = {0, W_star_nondeg_tm1[i]};
      //   Wstar_tm1_cens_temp = Wstar_tm1_cens;
      //   Wstar_tm1_cens_temp[i] = max(v);
      //   redist_nondeg = H*Wstar_tm1_cens_temp;
      //   Wt[i] = prop_deg*redist_deg[i] + (1-prop_deg)*redist_nondeg[i];
      // } else{
      //   Wt[i] = W_nondeg_t[i];
      // }
      Wt[i] = W_nondeg_t[i];
      W_star_t[i] = W_star_nondeg_t[i];
      w_tilde_prop = tnormRcpp(0, R_PosInf, Wt[i], w_prop_sd);
      lh1 = 0;
      lh3 = log(dtnorm(Wt[i], 0 ,  R_PosInf, w_tilde_prop, w_prop_sd))-
        log(dtnorm(w_tilde_prop, 0, R_PosInf, Wt[i], w_prop_sd));
      // lh3 = 0;
    }
    if(deg_temp[i] == 0){
      Wt[i] = W_nondeg_t[i];
      W_star_t[i] = W_star_nondeg_t[i];
      w_tilde_prop = tnormRcpp(0, R_PosInf, Wt[i], w_prop_sd);
      lh1 = R::dlnorm(Yt[i], w_tilde_prop, sqrt(s2_eta), TRUE) -
        R::dlnorm(Yt[i], Wt[i], sqrt(s2_eta), TRUE);
      lh3 = log(dtnorm(Wt[i], 0 ,  R_PosInf, w_tilde_prop, w_prop_sd))-
        log(dtnorm(w_tilde_prop, 0, R_PosInf, Wt[i], w_prop_sd));
    }
    W_tilde_prop = Wt;
    W_tilde_prop[i] = w_tilde_prop;
    
    for(int nn = 0; nn < ng; nn++){
      Ut.row(nn) = kron(Wt.row(nn).t(), Wt.row(nn).t(), U_ids);
      Vt.row(nn) = kron(Wt.row(nn).t(), XVt.col(nn), V_ids);
      Uprop_t.row(nn) = kron(W_tilde_prop.row(nn).t(), W_tilde_prop.row(nn).t(), U_ids);
      Vprop_t.row(nn) = kron(W_tilde_prop.row(nn).t(), XVt.col(nn), V_ids);
    }
    mean_curr = Wt + Vt * P + Ut * A;
    mean_prop = W_tilde_prop + Vprop_t * P + Uprop_t*A;
    for(int n = 0; n < ng*J; n++){
      lh2 = lh2+
        R::dnorm( W_star_t[n], mean_prop[n], s_vec[n], TRUE)-
        R::dnorm( W_star_t[n], mean_curr[n], s_vec[n], TRUE);
    } 
    
    if(deg_temp[i] == 1){
      lhr = 1000; // automatic accept
    } else{
      lhr = lh1 + lh2 + lh3;
    }
    
    u = runif(1)[0];
    if(w_tilde_prop >= w_lb & w_tilde_prop <= w_ub){
      if(lhr > log(u)){
        Wt[i] = w_tilde_prop;
        Ut = Uprop_t;
        Vt = Vprop_t;
        acc[i] = 1;
      }
    }
  }
  return List::create(Wt, deg_temp, Ut, Vt, acc);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_t_tilde3(const List& W_tilde, const List& W_star,  const List& W_star_deg,
                    const List& W_star_nondeg, const List& deg_id_ls,
                    const List& z_ids_ls, const List& Y_ls, const List& W_tilde_nondeg,
                    int ng, int J, double w_prop_sd, arma::vec s_vec,
                    double s2_eta, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                    arma::mat W0_tilde,arma::mat draw_deg_mat, int g,
                    const List& V_ls, const List& U_ls, const List& XV, 
                    double w_lb, double w_ub,
                    sp_mat H, arma::mat P, arma::mat A,
                    IntegerVector order, bool oneH = true) {
  IntegerVector ord = order;
  arma::mat Yt = Y_ls[t];
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat deg_temp = deg_id_ls[t];
  arma::vec z_ids = z_ids_ls[t];
  arma::mat W_nondeg_t = W_tilde_nondeg[t];
  int i, j;
  double w_tilde_prop;
  double lh1 = 0, lh2 = 0, lhr = 0, u; //, prop_deg;
  arma::vec acc = zeros(ng*J);
  arma::mat redist_deg, redist_nondeg, Wstar_tm1_cens_temp;
  NumericVector v(2);
  
  arma::mat Wstar_tm1_cens = W_star[t-1];
  arma::mat W_star_deg_tm1 = W_star_deg[t-1];
  arma::mat W_star_nondeg_tm1 = W_star_nondeg[t-1];
  arma::uvec neg_ids = find(Wstar_tm1_cens < 0);
  Wstar_tm1_cens.elem(neg_ids).fill(0);

  for(int q = 0; q < ng*J; q++){
    i = ord[q];
    
    j = ceil((i+1)/ng) - 1;
    
    if(deg_temp[i] == 1){
      w_tilde_prop = 0;
      lhr = 1000; //automatic accept
    }
    if(deg_temp[i] == 2){
      Wt[i] = W_nondeg_t[i];
      w_tilde_prop = tnormRcpp(0, R_PosInf, Wt[i], w_prop_sd);
      W_tilde_prop = Wt;
      lh1 = 0;
      lh2 = log(dtnorm(Wt[i], 0 ,  R_PosInf, w_tilde_prop, w_prop_sd))-
        log(dtnorm(w_tilde_prop, 0, R_PosInf, Wt[i], w_prop_sd));    
      lhr = lh1 + lh2;
    }
    if(deg_temp[i] == 0){
      Wt[i] =  W_nondeg_t[i];
      w_tilde_prop = tnormRcpp(0, R_PosInf, Wt[i], w_prop_sd);
      W_tilde_prop = Wt;
      W_tilde_prop[i] = w_tilde_prop;
      lh1 = R::dlnorm(Yt[i], w_tilde_prop, sqrt(s2_eta), TRUE) -
        R::dlnorm(Yt[i], Wt[i], sqrt(s2_eta), TRUE);
      lh2 = log(dtnorm(Wt[i], 0 ,  R_PosInf, w_tilde_prop, w_prop_sd))-
        log(dtnorm(w_tilde_prop, 0, R_PosInf, Wt[i], w_prop_sd));
      lhr = lh1 + lh2;
    }
    u = runif(1)[0];
    if(w_tilde_prop >= w_lb & w_tilde_prop <= w_ub){
      if(lhr > log(u) ){
        Wt[i] = w_tilde_prop;
        acc[i] = 1;
      }
    }
  }
  
  return List::create(Wt, deg_temp, acc);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_star_2level(const List& W_tilde, const List& W_star, const List& W_star_deg,
                const List& W_star_nondeg, const List& deg_id_ls, 
                const List& z_ids_ls, const List& Y_ls, arma::mat W_star_true_t,
                int ng, int J, arma::mat w_prop_sd, arma::vec s_vec,
                double s2_eta, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                arma::mat W0_tilde, 
                const List& V_ls, const List& U_ls, const List& XV,
                arma::mat wstar_ub, arma::mat wstar_lb,
                sp_mat H, arma::mat P, arma::mat A, 
                IntegerVector order, int ng2, arma::vec o, arma::uvec o_uvec, 
                const List& rel_eff_ls, bool oneH = true) {
  IntegerVector ord = order;
  arma::vec all_props = zeros(ng*J);
  arma::vec LH1 = zeros(ng*J), LH2 = zeros(ng*J), LH3= zeros(ng*J), LH4= zeros(ng*J);
  arma::mat Ytp1 = Y_ls[t+1];
  arma::mat W_tp1= W_tilde[t+1];
  arma::mat deg_tp1 = deg_id_ls[t+1];
  arma::uvec deg_id_tp1 = find( deg_tp1 == 1);
  arma::mat deg_t = deg_id_ls[t];
  arma::vec z_tp1_ids = z_ids_ls[t+1];
  arma::mat Utp1 = U_ls[t+1];
  arma::mat Vtp1 = V_ls[t+1];
  arma::mat Wstar_tp1 = W_star[t+1];
  
  arma::mat XVtp1 = XV[t+1];
  
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat Wstar_t = W_star[t];
  arma::mat Wstar_nondeg_t = W_star_nondeg[t];
  arma::mat Wstar_deg_t = W_star_deg[t];
  arma::mat deg_temp = deg_id_ls[t];
  arma::vec acc = zeros(ng*J);
  
  arma::mat Ut = U_ls[t];
  arma::mat Vt = V_ls[t];
  arma::mat XVt = XV[t];
  arma::mat m_t = Wt + Vt* P + Ut*A, mprop_tp1, m_tp1;
  int i, j;
  double lh1 = 0, lh2 = 0, lh3 = 0, lh4 = 0, lhr = 0, u;
  bool is_zero;
  arma::uvec neg_ids1, neg_ids2;
  double wstar_prop_t, wstar_t;
  arma::mat Wstar_prop_t, Wprop_tp1, Wstar_cens_prop, Uprop_tp1 = Utp1, Vprop_tp1 = Vtp1;
  arma::mat Wstar_cens, upscale_prop, upscale, redist_prop(ng2,J), redist(ng2,J);
  arma::mat redist_prop2(ng,J), redist2(ng, J);
  arma::vec rel_eff_tp1_vec= rel_eff_ls[t+1];
  arma::mat rel_eff_tp1(ng,J);
  for (int b = 0; b< J; b++){
    rel_eff_tp1.col(b) = rel_eff_tp1_vec;
  }
  // neg_ids2 = find(Wstar_cens < 0);
  // Wstar_cens.elem(neg_ids2).fill(0);
  // W_tp1 = H * Wstar_cens;
  // W_tp1.elem(deg_id_tp1).fill(0);
  // for(int nn = 0; nn < ng; nn++){
  //   Utp1.row(nn) = kron(W_tp1.row(nn).t(), W_tp1.row(nn).t(), U_ids);
  //   Vtp1.row(nn) = kron(W_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
  // }
  std::random_shuffle(ord.begin(), ord.end());
  for(int q = 0; q < ng*J; q++){
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    lh4 = 0;
    i = ord[q];
    
    // k = mod(i+1, ng);
    // if(k == 0){
    //   k = ng ;
    // }
    // k = k-1;
    j = ceil((i+1)/ng) - 1;
    
    if(deg_t[i] == 1){
      // when W_tilde = deg. 0
      wstar_t = Wstar_deg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = tnormRcpp(-0.3, 0.3, wstar_t, 0.075);
      lh4 = log(dtnorm(wstar_t, -0.3, 0.3, wstar_prop_t, 0.075))-
        log(dtnorm(wstar_prop_t, -0.3, 0.3, wstar_t, 0.075));
    } else{
      wstar_t = Wstar_nondeg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = R::rnorm(wstar_t, w_prop_sd[i]);
      lh4 = 0;
    }
    if(wstar_prop_t <= wstar_ub[i] & wstar_prop_t >= wstar_lb[i]){
      Wstar_prop_t = Wstar_t;
      Wstar_prop_t[i] = wstar_prop_t;
      // growth
      lh1 = R::dnorm(wstar_prop_t, m_t[i], s_vec[i], TRUE)-
        R::dnorm(wstar_t, m_t[i], s_vec[i], TRUE);
      
      // project forward
      Wstar_cens_prop = Wstar_prop_t;
      neg_ids1 = find(Wstar_cens_prop < 0);
      Wstar_cens_prop.elem(neg_ids1).fill(0);
      
      Wstar_cens = Wstar_t;
      neg_ids2 = find(Wstar_cens < 0);
      Wstar_cens.elem(neg_ids2).fill(0);
      
      // upscale
      upscale_prop = zeros(ng2,J);
      upscale = zeros(ng2,J);
      for(int r = 0; r < J; r++){
        upscale_prop.col(r) = tapply_sum_cpp(Wstar_cens_prop.col(r), o, ng2);
        upscale.col(r) = tapply_sum_cpp(Wstar_cens.col(r), o, ng2);
      }
      //redist: ng2 x J
      if(oneH){
        redist_prop = H*upscale_prop;
        redist = H * upscale;
      } else{
        redist_prop = H * vec_cpp(upscale_prop);
        redist = H * vec_cpp(upscale);
      }

      // redist2: ng x J
      for(int r = 0; r < J; r++){
        arma::mat temp_prop = redist_prop.col(r);
        arma::mat temp = redist.col(r);

        redist_prop2.col(r) = temp_prop.elem(o_uvec);
        redist2.col(r) = temp.elem(o_uvec);
      }

      // downscale
      Wprop_tp1 = redist_prop2 % rel_eff_tp1;
      W_tp1 = redist2 % rel_eff_tp1;

      // add in degenerate zeros
      Wprop_tp1.elem(deg_id_tp1).fill(0);
      W_tp1.elem(deg_id_tp1).fill(0);
      // 
      // Y_tp1 ~ LN(W_tile_tp1)
      for(int n = 0; n < ng*J; n++){
        is_zero = in_cpp(z_tp1_ids,n+1); 
        if(is_zero == false){
          lh2 = lh2 + R::dlnorm(Ytp1[n], Wprop_tp1[n], sqrt(s2_eta), TRUE)-
            R::dlnorm(Ytp1[n], W_tp1[n], sqrt(s2_eta), TRUE);
        }
      }
      
      for(int nn = 0; nn < ng; nn++){
        Utp1.row(nn) = kron(W_tp1.row(nn).t(), W_tp1.row(nn).t(), U_ids);
        Vtp1.row(nn) = kron(W_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
        Uprop_tp1.row(nn) = kron(Wprop_tp1.row(nn).t(), Wprop_tp1.row(nn).t(), U_ids);
        Vprop_tp1.row(nn) = kron(Wprop_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
      }
      
      mprop_tp1 =  Wprop_tp1 + Vprop_tp1 * P + Uprop_tp1 * A;
      m_tp1 =  W_tp1 + Vtp1* P + Utp1*A;
      for(int l = 0; l < ng*J; l++){
        lh3 = lh3 + R::dnorm(Wstar_tp1[l], mprop_tp1[l], s_vec[l], TRUE)-
          R::dnorm(Wstar_tp1[l], m_tp1[l], s_vec[l], TRUE);
      }
      
      lhr = lh1 +lh2 + lh3 + lh4;// + 
      // R::dnorm(wstar_prop_t,W_star_true_t[i], 0.05, TRUE)-
      // R::dnorm(wstar_t,W_star_true_t[i], 0.05, TRUE);
      u = runif(1)[0];
      if(lhr > log(u)){
        Wstar_t[i] = wstar_prop_t;
        W_tp1 = Wprop_tp1;
        Utp1 = Uprop_tp1;
        Vtp1 = Vprop_tp1;
        acc[i] = 1;
      }
    }
    LH1(i) = lh1;
    LH2(i) = lh2;
    LH3(i) = lh3;
    LH4(i) = lh4;
    all_props(i) = wstar_prop_t;
  }
  return List::create(Wstar_t, W_tp1, Utp1, Vtp1, acc, LH1, LH2, LH3, LH4, all_props,
                      m_t);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_TTm1_star_2level(const List& W_tilde, const List& W_star, const List& W_star_deg,
                     const List& W_star_nondeg, const List& deg_id_ls, 
                     const List& z_ids_ls, const List& Y_ls, arma::mat W_star_true_t,
                     int ng, int J, arma::mat w_prop_sd, arma::vec s_vec,
                     double s2_eta, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                     arma::mat W0_tilde,
                     const List& V_ls, const List& U_ls, const List& XV,
                     arma::mat wstar_ub, arma::mat wstar_lb,
                     sp_mat H, arma::mat P, arma::mat A, 
                     IntegerVector order, int ng2, arma::vec o, arma::uvec o_uvec, 
                     const List& rel_eff_ls, bool oneH = true) {
  IntegerVector ord = order;
  arma::mat Yt = Y_ls[t];
  arma::mat Ytp1 = Y_ls[t+1];
  arma::mat W_tp1 = W_tilde[t+1];
  arma::mat deg_tp1 = deg_id_ls[t+1];
  arma::uvec deg_id_tp1 = find( deg_tp1 == 1);
  arma::mat deg_t = deg_id_ls[t];
  arma::vec z_tp1_ids = z_ids_ls[t+1];
  
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat Wstar_t = W_star[t];
  arma::mat Wstar_deg_t = W_star_deg[t];
  arma::mat Wstar_nondeg_t = W_star_nondeg[t];
  arma::mat deg_temp = deg_id_ls[t];
  arma::vec acc = zeros(ng*J);
  // 
  arma::mat Ut = U_ls[t];
  arma::mat Vt = V_ls[t];
  arma::mat XVt = XV[t];
  arma::mat m_t = Wt + Vt* P + Ut*A;
  int i, j;
  arma::vec LH1 = zeros(ng*J), LH2 = zeros(ng*J), LH3= zeros(ng*J), LH4= zeros(ng*J),
    all_props = zeros(ng*J);
  double lh1 = 0, lh2 = 0, lh3 = 0, lh4 = 0, lhr = 0, u;
  bool is_zero;
  
  double wstar_prop_t, wstar_t;
  arma::mat Wstar_prop_t, Wprop_tp1, Wstar_cens_prop;
  
  arma::mat Wstar_cens, upscale_prop, upscale, redist_prop(ng2,J), redist(ng2,J);
  arma::mat redist_prop2(ng,J), redist2(ng, J);
  arma::vec rel_eff_tp1_vec= rel_eff_ls[t+1];
  arma::mat rel_eff_tp1(ng,J);
  for (int b = 0; b< J; b++){
    rel_eff_tp1.col(b) = rel_eff_tp1_vec;
  }
  arma::uvec neg_ids1, neg_ids2 ;
  std::random_shuffle(ord.begin(), ord.end());
  for(int q = 0; q < ng*J; q++){
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    lh4 = 0;
    i = ord[q];
    
    j = ceil((i+1)/ng) - 1;
    
    if(deg_t[i] == 1){
      // Wtilde = deg 0
      wstar_t = Wstar_deg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = tnormRcpp(-0.3, 0.3, wstar_t, 0.075);
      lh4 = log(dtnorm(wstar_t, -0.3, 0.3, wstar_prop_t, 0.075))-
        log(dtnorm(wstar_prop_t, -0.3, 0.3, wstar_t, 0.075));
    } else{
      wstar_t = Wstar_nondeg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = R::rnorm(wstar_t, w_prop_sd[i]);
    }
    if(wstar_prop_t <= wstar_ub[i] & wstar_prop_t >= wstar_lb[i]){
      Wstar_prop_t = Wstar_t;
      Wstar_prop_t[i] = wstar_prop_t;
      
      lh1 = R::dnorm(wstar_prop_t, m_t[i], s_vec[i], TRUE)-
        R::dnorm(wstar_t, m_t[i], s_vec[i], TRUE);

      // project forward
      Wstar_cens_prop = Wstar_prop_t;
      neg_ids1 = find(Wstar_cens_prop < 0);
      Wstar_cens_prop.elem(neg_ids1).fill(0);
      
      Wstar_cens = Wstar_t;
      neg_ids2 = find(Wstar_cens < 0);
      Wstar_cens.elem(neg_ids2).fill(0);
      
      // upscale
      upscale_prop = zeros(ng2,J);
      upscale = zeros(ng2,J);
      for(int r = 0; r < J; r++){
        upscale_prop.col(r) = tapply_sum_cpp(Wstar_cens_prop.col(r), o, ng2);
        upscale.col(r) = tapply_sum_cpp(Wstar_cens.col(r), o, ng2);
      }
      //redist: ng2 x J
      if(oneH){
        redist_prop = H*upscale_prop;
        redist = H * upscale;
      } else{
        redist_prop = H * vec_cpp(upscale_prop);
        redist = H * vec_cpp(upscale);
      }
      
      // redist2: ng x J
      for(int r = 0; r < J; r++){
        arma::mat temp_prop = redist_prop.col(r);
        arma::mat temp = redist.col(r);
        
        redist_prop2.col(r) = temp_prop.elem(o_uvec);
        redist2.col(r) = temp.elem(o_uvec);
      }
      
      // downscale
      Wprop_tp1 = redist_prop2 % rel_eff_tp1;
      W_tp1 = redist2 % rel_eff_tp1;
      
      // add in degenerate zeros
      Wprop_tp1.elem(deg_id_tp1).fill(0);
      W_tp1.elem(deg_id_tp1).fill(0);
      // 
      //
      for(int n = 0; n < ng*J; n++){
        is_zero = in_cpp(z_tp1_ids,n+1);
        if(is_zero == false){
          lh2 = lh2 + R::dlnorm(Ytp1[n], Wprop_tp1[n], sqrt(s2_eta), TRUE)-
            R::dlnorm(Ytp1[n], W_tp1[n], sqrt(s2_eta), TRUE);
        }
      }
      
      lhr = lh1 +lh2 + lh4 ;//+
      // R::dnorm(wstar_prop_t,W_star_true_t[i], 0.075, TRUE)-
      // R::dnorm(wstar_t,W_star_true_t[i], 0.075, TRUE);
      u = runif(1)[0];
      if(lhr > log(u)){
        Wstar_t[i] = wstar_prop_t;
        W_tp1 = Wprop_tp1;
        acc[i] = 1;
      } 
    }
    LH1(i) = lh1;
    LH2(i) = lh2;
    all_props(i) = wstar_prop_t;
  }
  return List::create(Wstar_t, W_tp1, acc, LH1, LH2, all_props);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_star_miss(const List& W_tilde, const List& W_star, const List& W_star_deg,
                const List& W_star_nondeg, const List& deg_id_ls, 
                const List& z_ids_ls, const List& Y_ls, 
                int ng_t, int ng_tp1, int J, arma::mat w_prop_sd, const List& s_vec_ls,
                const List& s_eta_vec_ls, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                const List& V_ls, const List& U_ls, const List& XV,
                arma::mat wstar_ub, arma::mat wstar_lb, const List& H_ls,
                arma::mat P, arma::mat A, 
                IntegerVector order, const List& obs_ids_ls,
                bool oneH = true) {
  IntegerVector ord = order;
  arma::vec all_props = zeros(ng_t*J);
  arma::vec s_vec_t = s_vec_ls[t];
  arma::vec s_vec_tp1 = s_vec_ls[t+1];
  arma::vec s_eta_vec_tp1 = s_eta_vec_ls[t+1];
  arma::vec LH1 = zeros(ng_t*J), LH2 = zeros(ng_t*J), LH3= zeros(ng_t*J), LH4= zeros(ng_t*J);
  arma::mat Ytp1 = Y_ls[t+1];
  // arma::mat prior_ind_t = prior_ind_ls[t];
  // arma::mat w_prior_t = w_prior_ls[t];
  sp_mat H_t = H_ls[t];
  arma::mat W_tp1= W_tilde[t+1];
  arma::mat deg_tp1 = deg_id_ls[t+1];
  arma::uvec deg_id_tp1 = find( deg_tp1 == 1);
  arma::mat deg_t = deg_id_ls[t];
  arma::vec z_tp1_ids = z_ids_ls[t+1];
  arma::mat Utp1 = U_ls[t+1];
  arma::mat Vtp1 = V_ls[t+1];
  arma::mat Wstar_tp1 = W_star[t+1];
  
  arma::mat XVtp1 = XV[t+1];
  
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat Wstar_t = W_star[t];
  arma::mat Wstar_nondeg_t = W_star_nondeg[t];
  arma::mat Wstar_deg_t = W_star_deg[t];
  arma::mat deg_temp = deg_id_ls[t];
  arma::vec acc = zeros(ng_t*J);
  
  arma::mat Ut = U_ls[t];
  arma::mat Vt = V_ls[t];
  arma::mat XVt = XV[t];
  arma::mat m_t = Wt + Vt* P + Ut*A, mprop_tp1, m_tp1;
  int i, j;
  double lh1 = 0, lh2 = 0, lh3 = 0, lh4 = 0, lh5 = 0, lhr = 0, u;
  bool is_zero;
  arma::uvec neg_ids1, neg_ids2;
  double wstar_prop_t, wstar_t;
  arma::mat Wstar_prop_t, Wprop_tp1, Wstar_cens_prop, Uprop_tp1 = Utp1, Vprop_tp1 = Vtp1;
  arma::mat Wstar_cens;

  std::random_shuffle(ord.begin(), ord.end());
  for(int q = 0; q < ng_t*J; q++){
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    lh4 = 0;
    lh5 = 0;
    i = ord[q];
    
    j = ceil((i+1)/ng_t) - 1;
    
    if(deg_t[i] == 1){
      // when W_tilde = deg. 0
      wstar_t = Wstar_deg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = tnormRcpp(-0.3, 0.3, wstar_t, 0.075);
      lh4 = log(dtnorm(wstar_t, -0.3, 0.3, wstar_prop_t, 0.075))-
        log(dtnorm(wstar_prop_t, -0.3, 0.3, wstar_t, 0.075));
    } else{
      wstar_t = Wstar_nondeg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = R::rnorm(wstar_t, w_prop_sd[i]);
      lh4 = 0;
    }
    if(wstar_prop_t <= wstar_ub[i] & wstar_prop_t >= wstar_lb[i]){
      Wstar_prop_t = Wstar_t;
      Wstar_prop_t[i] = wstar_prop_t;
      // growth
      lh1 = R::dnorm(wstar_prop_t, m_t[i], s_vec_t[i], TRUE)-
        R::dnorm(wstar_t, m_t[i], s_vec_t[i], TRUE);
      
      // project forward
      Wstar_cens_prop = Wstar_prop_t;
      neg_ids1 = find(Wstar_cens_prop < 0);
      Wstar_cens_prop.elem(neg_ids1).fill(0);
      if(oneH){
        Wprop_tp1 = H_t * Wstar_cens_prop;
      } else{
        Wprop_tp1 = vec2mat(H_t * vec_cpp(Wstar_cens_prop), ng_tp1, J);
      }
      Wprop_tp1.elem(deg_id_tp1).fill(0);
      
      Wstar_cens = Wstar_t;
      neg_ids2 = find(Wstar_cens < 0);
      Wstar_cens.elem(neg_ids2).fill(0);
      if(oneH){
        W_tp1 = H_t * Wstar_cens;
      } else{
        W_tp1 = vec2mat(H_t * vec_cpp(Wstar_cens), ng_tp1, J);
      }
      W_tp1.elem(deg_id_tp1).fill(0);
      // 
      // Y_tp1 ~ LN(W_tile_tp1)
      for(int n = 0; n < ng_tp1*J; n++){
        is_zero = in_cpp(z_tp1_ids,n+1); 
        if(is_zero == false){
          lh2 = lh2 + R::dlnorm(Ytp1[n], Wprop_tp1[n], s_eta_vec_tp1[n], TRUE)-
            R::dlnorm(Ytp1[n], W_tp1[n], s_eta_vec_tp1[n], TRUE);
        }
      }
      
      for(int nn = 0; nn < ng_tp1; nn++){
        Utp1.row(nn) = kron(W_tp1.row(nn).t(), W_tp1.row(nn).t(), U_ids);
        Vtp1.row(nn) = kron(W_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
        Uprop_tp1.row(nn) = kron(Wprop_tp1.row(nn).t(), Wprop_tp1.row(nn).t(), U_ids);
        Vprop_tp1.row(nn) = kron(Wprop_tp1.row(nn).t(), XVtp1.col(nn), V_ids);
      }
      
      mprop_tp1 =  Wprop_tp1 + Vprop_tp1 * P + Uprop_tp1 * A;
      m_tp1 =  W_tp1 + Vtp1* P + Utp1*A;
      for(int l = 0; l < ng_tp1*J; l++){
        lh3 = lh3 + R::dnorm(Wstar_tp1[l], mprop_tp1[l], s_vec_tp1[l], TRUE)-
          R::dnorm(Wstar_tp1[l], m_tp1[l], s_vec_tp1[l], TRUE);
      }
      // if(prior_ind_t[i] == 1){
      //   lh5 = R::dnorm4(wstar_prop_t, w_prior_t[i], 0.1, TRUE) -
      //     R::dnorm4(wstar_t , w_prior_t[i], 0.1, TRUE);
      // }
      
      lhr = lh1 +lh2 + lh3 + lh4 + lh5;// + 
      // R::dnorm(wstar_prop_t,W_star_true_t[i], 0.05, TRUE)-
      // R::dnorm(wstar_t,W_star_true_t[i], 0.05, TRUE);
      u = runif(1)[0];
      if(lhr > log(u)){
        Wstar_t[i] = wstar_prop_t;
        W_tp1 = Wprop_tp1;
        Utp1 = Uprop_tp1;
        Vtp1 = Vprop_tp1;
        acc[i] = 1;
      }
    }
    LH1(i) = lh1;
    LH2(i) = lh2;
    LH3(i) = lh3;
    LH4(i) = lh4;
    all_props(i) = wstar_prop_t;
  }
  return List::create(Wstar_t, W_tp1, Utp1, Vtp1, acc, LH1, LH2, LH3, LH4, all_props,
                      m_t);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List sampW_TTm1_star_miss(const List& W_tilde, const List& W_star, const List& W_star_deg,
                     const List& W_star_nondeg, const List& deg_id_ls, 
                     const List& z_ids_ls, const List& Y_ls, 
                     int ng_t, int ng_tp1, int J, arma::mat w_prop_sd, const List& s_vec_ls,
                     const List& s_eta_vec_ls, int t, int TT, arma::uvec U_ids, arma::uvec V_ids,
                     const List& V_ls, const List& U_ls, const List& XV,
                     arma::mat wstar_ub, arma::mat wstar_lb, const List& H_ls,
                     arma::mat P, arma::mat A, 
                     IntegerVector order, bool oneH = true) {
  IntegerVector ord = order;
  arma::vec s_vec_t = s_vec_ls[t];
  arma::vec s_vec_tp1 = s_vec_ls[t+1];
  arma::vec s_eta_vec_tp1 = s_eta_vec_ls[t+1];
  sp_mat H = H_ls[t];
  arma::mat Yt = Y_ls[t];
  arma::mat Ytp1 = Y_ls[t+1];
  arma::mat W_tp1 = W_tilde[t+1];
  arma::mat deg_tp1 = deg_id_ls[t+1];
  arma::uvec deg_id_tp1 = find( deg_tp1 == 1);
  arma::mat deg_t = deg_id_ls[t];
  arma::vec z_tp1_ids = z_ids_ls[t+1];
  
  arma::mat Wt = W_tilde[t], W_tilde_prop;
  arma::mat Wstar_t = W_star[t];
  arma::mat Wstar_deg_t = W_star_deg[t];
  arma::mat Wstar_nondeg_t = W_star_nondeg[t];
  arma::mat deg_temp = deg_id_ls[t];
  arma::vec acc = zeros(ng_t*J);
  // 
  arma::mat Ut = U_ls[t];
  arma::mat Vt = V_ls[t];
  arma::mat XVt = XV[t];
  arma::mat m_t = Wt + Vt* P + Ut*A;
  int i, j;
  arma::vec LH1 = zeros(ng_t*J), LH2 = zeros(ng_t*J), LH3= zeros(ng_t*J), LH4= zeros(ng_t*J),
    all_props = zeros(ng_t*J);
  double lh1 = 0, lh2 = 0, lh3 = 0, lh4 = 0, lhr = 0, u;
  bool is_zero;
  
  double wstar_prop_t, wstar_t;
  arma::mat Wstar_prop_t, Wprop_tp1, Wstar_cens_prop;
  
  arma::mat Wstar_cens;
  arma::uvec neg_ids1, neg_ids2 ;
  std::random_shuffle(ord.begin(), ord.end());
  for(int q = 0; q < ng_t*J; q++){
    lh1 = 0;
    lh2 = 0;
    lh3 = 0;
    lh4 = 0;
    i = ord[q];
    
    j = ceil((i+1)/ng_t) - 1;
    
    if(deg_t[i] == 1){
      // Wtilde = deg 0
      wstar_t = Wstar_deg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = tnormRcpp(-0.3, 0.3, wstar_t, 0.075);
      lh4 = log(dtnorm(wstar_t, -0.3, 0.3, wstar_prop_t, 0.075))-
        log(dtnorm(wstar_prop_t, -0.3, 0.3, wstar_t, 0.075));
    } else{
      wstar_t = Wstar_nondeg_t[i];
      Wstar_t[i] = wstar_t;
      wstar_prop_t = R::rnorm(wstar_t, w_prop_sd[i]);
    }
    if(wstar_prop_t <= wstar_ub[i] & wstar_prop_t >= wstar_lb[i]){
      Wstar_prop_t = Wstar_t;
      Wstar_prop_t[i] = wstar_prop_t;
      
      lh1 = R::dnorm(wstar_prop_t, m_t[i], s_vec_t[i], TRUE)-
        R::dnorm(wstar_t, m_t[i], s_vec_t[i], TRUE);
      Wstar_cens_prop = Wstar_prop_t;
      neg_ids1 = find(Wstar_cens_prop < 0);
      Wstar_cens_prop.elem(neg_ids1).fill(0);
      if(oneH){
        Wprop_tp1 = H * Wstar_cens_prop;
      } else{
        Wprop_tp1 = vec2mat(H * vec_cpp(Wstar_cens_prop),ng_tp1,J);
      }
      Wprop_tp1.elem(deg_id_tp1).fill(0);
      
      Wstar_cens = Wstar_t;
      neg_ids2 = find(Wstar_cens < 0);
      Wstar_cens.elem(neg_ids2).fill(0);
      if(oneH){
        W_tp1 = H * Wstar_cens;
      } else{
        W_tp1 = vec2mat(H * vec_cpp(Wstar_cens),ng_tp1, J);
      }
      W_tp1.elem(deg_id_tp1).fill(0);
      //
      for(int n = 0; n < ng_tp1*J; n++){
        is_zero = in_cpp(z_tp1_ids,n+1);
        if(is_zero == false){
          lh2 = lh2 + R::dlnorm(Ytp1[n], Wprop_tp1[n], s_eta_vec_tp1[n], TRUE)-
            R::dlnorm(Ytp1[n], W_tp1[n], s_eta_vec_tp1[n], TRUE);
        }
      }
      
      lhr = lh1 +lh2 + lh4 ;//+
      // R::dnorm(wstar_prop_t,W_star_true_t[i], 0.075, TRUE)-
      // R::dnorm(wstar_t,W_star_true_t[i], 0.075, TRUE);
      u = runif(1)[0];
      if(lhr > log(u)){
        Wstar_t[i] = wstar_prop_t;
        W_tp1 = Wprop_tp1;
        acc[i] = 1;
      } 
    }
    LH1(i) = lh1;
    LH2(i) = lh2;
    all_props(i) = wstar_prop_t;
  }
  return List::create(Wstar_t, W_tp1, acc, LH1, LH2, all_props);
}





