functions {
  // anisotropic squared exponential covariance
  matrix cov_sqExp_anis(matrix x,
                        real alpha,
                        vector rho,
                        real delta) {
    int N = rows(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        // K[i, j] = sq_alpha* exp(-0.5 * dot_self((x[i,] - x[j,]) ./ rho));
        K[i, j] = sq_alpha* exp(-0.5 * dot_self((to_vector(x[i,]) - to_vector(x[j,])) ./ rho));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return K;
  }
  // covariance between computer model data and field data
  matrix cov_off_diag(matrix x1,
                      matrix x2,
                      real alpha,
                      vector rho) {
    int N1 = rows(x1);
    int N2 = rows(x2);
    matrix[N1, N2] K; 
    real sq_alpha = square(alpha);
    for (i in 1:N1) {
      for (j in 1:N2) {
        K[i, j] = sq_alpha* exp(-0.5 * dot_self((to_vector(x1[i,] )- to_vector(x2[j,])) ./ rho));
      }
    }
    return K;
  }
  // Kennedy and O'Hagan (2001) covariance matrix
    matrix KOH_cov(matrix xt_M,
                 matrix xt_F,
                 matrix x_F,
                 real alpha_M,
                 vector rho_M,
                 real alpha_B,
                 vector rho_B,
                 real sigma,
                 real delta){
    int m = rows(xt_M);
    int n = rows(xt_F);
    matrix[n+m, n+m] K;
    // matrix[m, n] K_off_diag;
    // K_off_diag = cov_off_diag(xt_M, xt_F, alpha_M, rho_M);
    matrix[n, m] K_off_diag;
    K_off_diag = cov_off_diag(xt_F, xt_M, alpha_M, rho_M);
    K[1:n,1:n] = cov_sqExp_anis(xt_F, alpha_M, rho_M, delta) + cov_sqExp_anis(x_F, alpha_B, rho_B, square(sigma));// K[1:m,1:m] = cov_sqExp_anis(xt_M, alpha_M, rho_M,delta);
    K[1:n, (n+1):(m+n)] = K_off_diag;
    K[(n+1):(m+n), 1:n] = K_off_diag';
    K[(n+1):(m+n), (n+1):(m+n)] = cov_sqExp_anis(xt_M, alpha_M, rho_M,delta);
    
    return cholesky_decompose(K);
  }
  // add disscrepancy in the predictionss
  vector Pgp_pred_rng(matrix xt_M,
                      vector t_F,
                      matrix x_F,
                      matrix x_Fpred,
                      vector y,
                      vector eta,
                      real alpha_M,
                      vector rho_M,
                      real alpha_B,
                      vector rho_B,
                      real sigma,
                      real delta) {

    int m = rows(xt_M);
    int n = rows(x_F);
    int p = cols(x_F);
    int q = cols(xt_M)-cols(x_F);
    int N = m + n;
    int n_pred = rows(x_Fpred);
    vector[n_pred] f2;
    {
      matrix[N, N] L_K;
      vector[N] K_div_y;
      matrix[n_pred, N] qT_p;
      matrix[N, n_pred] v_pred;
      vector[n_pred] f2_mu;
      matrix[n_pred, n_pred] cov_f2;
      matrix[n_pred, n_pred] diag_delta;
      
      vector[N] y_eta;
      matrix[n, p+q] xt_F;
      matrix[n_pred, p+q] xt_Fpred;
      
      y_eta[1:n] = y;
      y_eta[(n+1):N] = eta;
      xt_F = append_col(x_F, rep_matrix(to_row_vector(t_F), n));
      xt_Fpred = append_col(x_Fpred, rep_matrix(to_row_vector(t_F), n_pred));
      
      L_K =  KOH_cov(xt_M, xt_F, x_F, alpha_M, rho_M, alpha_B, rho_B, sigma, delta);// L_K =  K_wk2(tP, tI, rho, alpha, sigmaP, sigmaI, R, C);
      K_div_y = mdivide_left_tri_low(L_K, y_eta);
      K_div_y = mdivide_right_tri_low(K_div_y', L_K)';
      // continue here
      // remember to create xt_Fpred, from x_Fpred and t_F
      qT_p[1:n_pred, 1:n] = cov_off_diag(xt_Fpred, xt_F, alpha_M, rho_M)+ cov_off_diag(x_Fpred, x_F, alpha_B, rho_B);// qT_p[1:n_pred, 1:n] = KP_pred(tP, x_Fpred, rho, alpha);
      qT_p[1:n_pred, (n+1):N] = cov_off_diag(xt_Fpred, xt_M, alpha_M, rho_M);// qT_p[1:n_pred, (nP+1):N] = KPI_pred(tI, x_Fpred, rho, alpha, R, C);
      f2_mu = (qT_p * K_div_y);
      v_pred = mdivide_left_tri_low(L_K, qT_p');
      cov_f2 = cov_sqExp_anis(xt_Fpred, alpha_M, rho_M, delta) + cov_sqExp_anis(x_Fpred, alpha_B, rho_B, square(sigma))- v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(1e-6, n_pred));
                                                                      
      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
}

data {
  // define computer model and field data
  int<lower=1> m; // # of computer model data
  int<lower=1> n; // # of field data
  int<lower=1> p; // # of observable inputs 
  int<lower=1> q; // # of calibration parameters
  
  vector[m] eta;   // computer simulation outputs
  vector[n] y;     // field outputs
  
  // vector[D] x[N];
  matrix[m,p] x_M; // computer model design inputs
  matrix[m,q] t_M; // computer model design parameters
  matrix[n,p] x_F; // field inputs
  real delta; //nugget
  
  int<lower=1> n_pred;
  matrix[n_pred,p] x_Fpred;
  int<lower=1> N_samples;
  matrix[N_samples,p+q] rho_M;
  matrix[N_samples,p] rho_B;
  vector[N_samples] alpha_M;
  vector[N_samples] alpha_B;
  vector[N_samples] sigma;
  
  matrix[N_samples,q] t_F;
}

transformed data {
  vector[m+n] y_eta;
  matrix[m, p+q] xt_M; 
  y_eta = append_row(y, eta);
  xt_M = append_col(x_M, t_M);
}

parameters {}
model {}
generated quantities {
  vector[n_pred] f_P;
  matrix[N_samples, n_pred] y_P;

  for(ns in 1:N_samples) {
    //f_P = Pgp_pred_rng(tP, tI, yP, yI, x_Fpred[1:n_pred], alpha[n], rho[n], rho_d[n], alpha_d[n], R[n], C[n], sigmaP[n], sigmaI[n]);
    f_P = Pgp_pred_rng(xt_M, to_vector(t_F[ns,]), x_F, x_Fpred, y, eta, alpha_M[ns], to_vector(rho_M[ns,]), alpha_B[ns], to_vector(rho_B[ns,]), sigma[ns], delta);
    for(np in 1:n_pred) {
      y_P[ns, np] = normal_rng(f_P[np], sigma[ns]);
    }  
  }
  
}
