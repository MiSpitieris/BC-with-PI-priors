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
    matrix[n, m] K_off_diag;
    K_off_diag = cov_off_diag(xt_F, xt_M, alpha_M, rho_M);
    K[1:n,1:n] = cov_sqExp_anis(xt_F, alpha_M, rho_M, delta) + cov_sqExp_anis(x_F, alpha_B, rho_B, square(sigma));
    K[1:n, (n+1):(m+n)] = K_off_diag;
    K[(n+1):(m+n), 1:n] = K_off_diag';
    K[(n+1):(m+n), (n+1):(m+n)] = cov_sqExp_anis(xt_M, alpha_M, rho_M,delta);
    
    return cholesky_decompose(K);
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
}

transformed data {
  vector[m+n] y_eta; 
  y_eta = append_row(y, eta);
}

parameters {
  real<lower=0> alpha_M;
  real<lower=0> alpha_B;
  vector<lower=0>[p+q] rho_M;
  vector<lower=0>[p] rho_B;
  real<lower=0,upper=10> sigma;
  real mu_M;
  row_vector<lower=0.5, upper=3>[q] t_F;
}

model {
  matrix[m+n, m+n] L_K;
  matrix[m, p+q] xt_M; 
  matrix[n, p+q] xt_F;
  vector[m+n] mu; 

  mu[1:n] = rep_vector(mu_M, n); 
  mu[(n+1): (m+n)] = rep_vector(mu_M, m);
  
  xt_M = append_col(x_M, t_M);
  xt_F = append_col(x_F, rep_matrix(t_F, n));
  
  rho_M ~ gamma(3.0/2, 1.0);
  rho_B ~ gamma(3.0/2, 1.0);
  // alpha_M ~ inv_gamma(5, 5);
  // alpha_B ~ inv_gamma(5, 5);
  alpha_M ~ normal(sd(eta),4);
  alpha_B ~ normal(10,50);
  mu_M ~ normal(mean(eta),10);
  
  L_K = KOH_cov(xt_M, xt_F, x_F, alpha_M, rho_M, alpha_B, rho_B, sigma, delta);
  y_eta ~ multi_normal_cholesky(mu, L_K);
}

