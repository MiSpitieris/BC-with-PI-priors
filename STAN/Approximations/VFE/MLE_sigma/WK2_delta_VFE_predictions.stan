// predictions for FITC models
functions {
  // Physics informed prior mean of the WK2 model
  vector mu_fn(real mu_wk2, 
               real R, 
               int nP, 
               int nI){
                 vector[nP]  mP = rep_vector(mu_wk2,nP); 
                 vector[nI]  mI = rep_vector((1/R)*mu_wk2,nI);
                 vector[nP + nI] mu;
                 mu = append_row(mP, mI);
                 return(mu);
  }
  
  matrix K_PP(vector tP_i,
              vector tP_j,
              real rho,
              real alpha,
              real rho_d,
              real alpha_d){
    int nP_i = rows(tP_i);
    int nP_j = rows(tP_j);
    matrix[nP_i, nP_j] K;
    matrix[nP_i, nP_j] KB;
    real alpha2 = pow(alpha, 2.0);
    real rho2_ = pow(rho, -2.0);
    
    for (i in 1:nP_i){
      for (j in 1:nP_j){
        K[i,j] = exp(-pow(tP_i[i] - tP_j[j], 0.2e1) * rho2_);
        K[i,j] = alpha2 * K[i,j];
      }
    }
    
    for (i in 1:nP_i){
      for (j in 1:nP_j){
        KB[i,j] = exp(-pow(tP_i[i] - tP_j[j], 0.2e1) * pow(rho_d, -0.2e1));
        KB[i,j] = pow(alpha_d, 0.2e1) * KB[i,j];
      }
    }
    K[1:nP_i, 1:nP_j] = K[1:nP_i, 1:nP_j] + KB[1:nP_i, 1:nP_j];
    return(K);
   }
   
   matrix K_PI(vector tP, 
               vector tI,
               real rho,
               real alpha,
               real R,
               real C){
    int nP = rows(tP);
    int nI = rows(tI);
    matrix[nP, nI] K;
    real alpha2 = pow(alpha, 2.0);
    real rho2_ = pow(rho, -2.0);
    
    for (i in 1:nP){
      for (j in 1:nI){
        K[i, j] = 0.1e1 / R * exp(-pow(tP[i] - tI[j], 0.2e1) * rho2_) 
        + 0.2e1 * C * (tP[i] - tI[j])* rho2_ 
        * exp(-pow(tP[i] - tI[j], 0.2e1) * rho2_);
       K[i, j] = alpha2 * K[i,j];              
      }
   }
   return(K);
  }
  
   matrix K_IP(vector tI, 
               vector tP,
               real rho,
               real alpha,
               real R,
               real C){
    int nP = rows(tP);
    int nI = rows(tI);
    matrix[nI, nP] K;
    real alpha2 = pow(alpha, 2.0);
    real rho2_ = pow(rho, -2.0);
    
        for (i in 1:nI){
          for (j in 1:nP){
            K[i, j] = 0.1e1 / R * exp(-pow(tI[i] - tP[j], 0.2e1) * rho2_)
            - 0.2e1 * C * (tI[i] - tP[j]) * rho2_
            * exp(-pow(tI[i] - tP[j], 0.2e1) * rho2_);
            K[i, j] = alpha2 * K[i, j];
          }
        }
        return(K);
  }
  
  matrix K_II(vector tI_i, 
              vector tI_j,
              real rho,
              real alpha,
              real R,
              real C){
    int nI_i = rows(tI_i);
    int nI_j = rows(tI_j);
    matrix[nI_i, nI_j] K;
    real alpha2 = pow(alpha, 2.0);
    real rho2_ = pow(rho, -2.0);
    real rho4_ = pow(rho2_, 2.0);
    real R2_ = pow(R, -2.0);
    real C2 = pow(C, 2.0);
    
    for (i in 1:nI_i){
          for (j in 1:nI_j){
            K[ i, j] = 
              R2_ * exp(-pow(tI_i[i] - tI_j[j], 0.2e1) * rho2_) 
            + C2 * (0.2e1 * rho2_ * exp(-pow(tI_i[i] - tI_j[j], 0.2e1) 
                                                          * rho2_) - 0.4e1 * pow(tI_i[i] - tI_j[j], 0.2e1) 
                           * rho4_ * exp(-pow(tI_i[i] - tI_j[j], 0.2e1) * rho2_));
            K[i, j] = alpha2 * K[i, j];
          }
    }
  return(K);
  }
  
  matrix K_wk2_XZ(vector tP,
                  vector tI,
                  vector zP,
                  vector zI,
                  real rho,
                  real alpha,
                  real rho_d,
                  real alpha_d,
                  real R,
                  real C){
    int nP = rows(tP);
    int nI = rows(tI);
    int mP = rows(zP);
    int mI = rows(zI);
    matrix[nP+nI, mP+mI] K;
    K[1:nP, 1:mP] = K_PP(tP,zP, rho,alpha, rho_d, alpha_d); 
    K[1:nP, (mP+1):(mP+mI)] = K_PI(tP, zI, rho, alpha, R, C);
    K[(nP+1):(nP+nI), 1:mP] = K_IP(tI, zP, rho, alpha, R, C);
    K[(nP+1):(nP+nI),(mP+1):(mP+mI)] = K_II(tI, zI, rho, alpha, R, C);
    return(K);
  }  
  
  // Physics informed prior kernel of the WK2 model
  matrix K_wk2(vector tP,
               vector tI,
               real rho,
               real alpha,
               real rho_d,
               real alpha_d,
               real R,
               real C) {
    int nP = rows(tP);
    int nI = rows(tI);
    matrix[nP + nI, nP + nI] K;
    matrix[nP, nP] KB;

    real alpha2 = pow(alpha, 2.0);
    real rho2_ = pow(rho, -2.0);
    real rho4_ = pow(rho2_, 2.0);
    real R2_ = pow(R, -2.0);
    real C2 = pow(C, 2.0);
    // KP
for (i in 1:(nP-1)){
  K[i,i] = alpha2;
  for (j in (i+1):nP){
    K[i,j] = exp(-pow(tP[i] - tP[j], 0.2e1) * rho2_);
    K[i,j] = alpha2 * K[i,j];
    K[j,i] = K[i,j];
  }
  K[nP,nP] = alpha2;
}
    // K_delta (press bias)
    for (i in 1:(nP-1)){
      KB[i,i] = pow(alpha_d, 0.2e1);
      for (j in (i+1):nP){
        KB[i,j] = exp(-pow(tP[i] - tP[j], 0.2e1) * pow(rho_d, -0.2e1));
        KB[i,j] = pow(alpha_d, 0.2e1) * KB[i,j];
        KB[j,i] = KB[i,j];
      }
      KB[nP,nP] = pow(alpha_d, 0.2e1);
    }
    K[1:nP, 1:nP] = K[1:nP, 1:nP] + KB[1:nP, 1:nP];

   // KPI
   for (i in 1:nP){
     for (j in 1:nI){
       K[i, nP + j] = 0.1e1 / R * exp(-pow(tP[i] - tI[j], 0.2e1) * rho2_)
       + 0.2e1 * C * (tP[i] - tI[j])* rho2_
       * exp(-pow(tP[i] - tI[j], 0.2e1) * rho2_);
       K[i, nP + j] = alpha2 * K[i, nP + j];
     }
   }

   // KIP (KIP = KPI')
   K[(nP + 1):(nP + nI), 1:nP] = K[1:nP, (nP + 1):(nP + nI)]';

   // KI
   for (i in 1:(nI-1)){
     K[nP + i, nP +i] =
     R2_ * exp(-pow(tI[i] - tI[i], 0.2e1) * rho2_)
     + C2 * (0.2e1 * rho2_ * exp(-pow(tI[i] - tI[i], 0.2e1)* rho2_) 
     - 0.4e1 * pow(tI[i] - tI[i], 0.2e1)* rho4_ * exp(-pow(tI[i] - tI[i], 0.2e1) * rho2_));
     K[nP + i, nP +i] = alpha2 * K[nP + i, nP +i];
     for (j in (i+1):nI){
       K[nP + i, nP +j] =
       R2_ * exp(-pow(tI[i] - tI[j], 0.2e1) * rho2_)
       + C2 * (0.2e1 * rho2_ * exp(-pow(tI[i] - tI[j], 0.2e1)* rho2_) 
       - 0.4e1 * pow(tI[i] - tI[j], 0.2e1) * rho4_ * exp(-pow(tI[i] - tI[j], 0.2e1) * rho2_));
       K[nP + i, nP +j] = alpha2 * K[nP + i, nP +j];
       K[nP + j, nP +i] = K[nP + i, nP +j];
     }
     K[nP + nI, nP +nI] =
      R2_ * exp(-pow(tI[nI] - tI[nI], 0.2e1) * rho2_)
      + C2 * (0.2e1 * rho2_ * exp(-pow(tI[nI] - tI[nI], 0.2e1)* rho2_) 
      - 0.4e1 * pow(tI[nI] - tI[nI], 0.2e1)* rho4_ * exp(-pow(tI[nI] - tI[nI], 0.2e1) * rho2_));
      K[nP + nI, nP +nI] = alpha2 * K[nP + nI, nP +nI];
    }
    return K;
  }
  
  
  vector K_diag(vector tP,
               vector tI,
               real rho,
               real alpha,
               real rho_d,
               real alpha_d,
               real R,
               real C) {
    int nP = rows(tP);
    int nI = rows(tI);
    vector[nP + nI] v;

    real alpha2 = pow(alpha, 2.0);
    real alpha_d2 = pow(alpha_d, 2.0);
    real rho2_ = pow(rho, -2.0);
    real rho4_ = pow(rho2_, 2.0);
    real R2_ = pow(R, -2.0);
    real C2 = pow(C, 2.0);
    // KP
    for (i in 1:nP){
      v[i] = alpha2+alpha_d2;
    }

   for (i in 1:nI){
     v[nP + i] =
     R2_ * exp(-pow(tI[i] - tI[i], 0.2e1) * rho2_)
     + C2 * (0.2e1 * rho2_ * exp(-pow(tI[i] - tI[i], 0.2e1)* rho2_) 
     - 0.4e1 * pow(tI[i] - tI[i], 0.2e1)* rho4_ * exp(-pow(tI[i] - tI[i], 0.2e1) * rho2_));
     v[nP + i] = alpha2 * v[nP + i];
  }
  // K = diag_matrix(v);
  return v;
  }
  // pred P
  vector Pgp_pred_rng(vector tP,
                      vector tI,
                      vector yP,
                      vector yI,
                      vector tP_pred,
                      real alpha,
                      real rho,
                      real rho_d,
                      real alpha_d,
                      real R,
                      real C,
                      real sigmaP,
                      real sigmaI,
                      vector zP,
                      vector zI)  {
    int nP = rows(tP);
    int nI = rows(tI);
    int N = nP + nI;
    int mP = rows(zP);
    int mI = rows(zI);
    int m = mP +mI;
    int n_pred = rows(tP_pred);
    int s=n_pred;
    vector[n_pred] f_pred;
    {
      vector[N] y;
      matrix[N,m] K_fu;
      matrix[m,N] K_uf;
      // matrix[N,N] K_ff;
      vector[N] K_ff;
      matrix[m,m] K_uu;
      matrix[m,m] L_uu;
      matrix[m,n_pred] K_us;
      matrix[n_pred, m] K_su;
      matrix[n_pred, n_pred] K_ss;
      matrix[m,N] A;
      matrix[m,m] B;
      matrix[m,m] AS_invA;
      matrix[m,m] B_chol;
      matrix[m,1] v_mat;
      vector[n_pred] mu;
      vector[N] sigma_vec_inv;
      matrix[N,N] Sigma_inv;
      matrix[N,N] sigmaI_inv;
      vector[N] r;
      vector[N] r_l;
      vector[m] c;
      
      matrix[m,s] CC;
      matrix[m,s] D;
      matrix[s,s] cov;
      
      y[1:nP] = yP[1:nP];
      y[(nP+1):N] = yI[1:nI];
      
      K_fu = K_wk2_XZ(tP, tI, zP, zI, rho, alpha, rho_d, alpha_d, R, C);
      K_uf = K_fu';
      K_uu = K_wk2(zP, zI, rho, alpha, rho_d, alpha_d, R, C)+ diag_matrix(rep_vector(1e-6,m));
      K_ff = K_diag(tP, tI, rho, alpha, rho_d, alpha_d, R, C);
      K_ss = K_PP(tP_pred, tP_pred, rho, alpha, rho_d, alpha_d);
      K_su[1:n_pred,1:mP] =  K_PP(tP_pred, zP, rho, alpha, rho_d, alpha_d);
      K_su[1:n_pred,(mP+1):m] = K_PI(tP_pred, zI, rho, alpha, R, C);
      K_us = K_su';
      
      sigma_vec_inv = to_vector(append_row(rep_vector(1/square(sigmaP),nP),rep_vector(1/square(sigmaI),nI)));
      Sigma_inv = diag_matrix(sigma_vec_inv);
  
      L_uu = cholesky_decompose(K_uu);
      A = mdivide_left_tri_low(L_uu, K_uf);
      AS_invA = A*Sigma_inv*A';
      B = diag_matrix(rep_vector(1,m)) + AS_invA;
      B_chol = cholesky_decompose(B+ diag_matrix(rep_vector(1e-6,m)));
      
      // sigma_vec = to_vector(append_row(rep_vector(square(sigmaP),nP),rep_vector(square(sigmaI),nI)));
      // L_uu = cholesky_decompose(K_uu);
      // A = mdivide_left_tri_low(L_uu, K_uf);
      // Lambda = diag_matrix(K_ff-diagonal(A'*A) +sigma_vec);
      // Lambda_inv = Lambda;
      // for (i in 1:N){
      //   Lambda_inv[i,i] = 1.0/Lambda[i,i];
      // }
      // B = diag_matrix(rep_vector(1,m)) + A*Lambda_inv*A';
      // B_chol = cholesky_decompose(B+ diag_matrix(rep_vector(1e-6,m)));
      
      // r = y;
      // r_l = Lambda_inv*r;
      
      r = y;
      r_l = Sigma_inv*r;
      c = mdivide_left_tri_low(B_chol, A*r_l);
      c = mdivide_right_tri_low(c', B_chol)';
      c = mdivide_right_tri_low(c', L_uu)';
      mu = K_su* c;

      CC = mdivide_left_tri_low(L_uu, K_us);
      D = mdivide_left_tri_low(B_chol, CC);
      cov = K_ss - CC'*CC + D'*D + diag_matrix(rep_vector(1e-6, n_pred));
      f_pred = multi_normal_rng(mu, cov);
    }
    return f_pred;
  }
  
  vector Igp_pred_rng(vector tP,
                      vector tI,
                      vector yP,
                      vector yI,
                      vector tI_pred,
                      real alpha,
                      real rho,
                      real rho_d,
                      real alpha_d,
                      real R,
                      real C,
                      real sigmaP,
                      real sigmaI,
                      vector zP,
                      vector zI)  {
    int nP = rows(tP);
    int nI = rows(tI);
    int N = nP + nI;
    int mP = rows(zP);
    int mI = rows(zI);
    int m = mP +mI;
    int n_pred = rows(tI_pred);
    int s=n_pred;
    vector[n_pred] f_pred;
    {
      vector[N] y;
      matrix[N,m] K_fu;
      matrix[m,N] K_uf;
      vector[N] K_ff;
      matrix[m,m] K_uu;
      matrix[m,m] L_uu;
      matrix[m,n_pred] K_us;
      matrix[n_pred, m] K_su;
      matrix[n_pred, n_pred] K_ss;
      matrix[m,N] A;
      matrix[m,m] B;
      matrix[m,m] B_chol;
      matrix[m,1] v_mat;
      vector[n_pred] mu;
      vector[N] sigma_vec_inv;
      matrix[m,m] AS_invA;
      matrix[N,N] Sigma_inv;
      vector[N] r;
      vector[N] r_l;
      vector[m] c;

      matrix[m,s] CC;
      matrix[m,s] D;
      matrix[s,s] cov;
      
      y[1:nP] = yP[1:nP];
      y[(nP+1):N] = yI[1:nI];
      
      K_fu = K_wk2_XZ(tP, tI, zP, zI, rho, alpha, rho_d, alpha_d, R, C);
      K_uf = K_fu';
      K_uu = K_wk2(zP, zI, rho, alpha, rho_d, alpha_d, R, C)+ diag_matrix(rep_vector(1e-6,m));
      K_ff = K_diag(tP, tI, rho, alpha, rho_d, alpha_d, R, C);
      K_ss = K_II(tI_pred, tI_pred, rho, alpha, R, C);
      K_su[1:n_pred,1:mP] =  K_IP(tI_pred, zP, rho, alpha, R, C);
      K_su[1:n_pred,(mP+1):m] = K_II(tI_pred, zI, rho, alpha, R, C);
      K_us = K_su';
      
      sigma_vec_inv = to_vector(append_row(rep_vector(1/square(sigmaP),nP),rep_vector(1/square(sigmaI),nI)));
      Sigma_inv = diag_matrix(sigma_vec_inv);
  
      L_uu = cholesky_decompose(K_uu);
      A = mdivide_left_tri_low(L_uu, K_uf);
      AS_invA = A*Sigma_inv*A';
      B = diag_matrix(rep_vector(1,m)) + AS_invA;
      B_chol = cholesky_decompose(B+ diag_matrix(rep_vector(1e-6,m)));
      r = y;
      r_l = Sigma_inv*r;

      c = mdivide_left_tri_low(B_chol, A*r_l);
      c = mdivide_right_tri_low(c', B_chol)';
      c = mdivide_right_tri_low(c', L_uu)';

      mu = K_su* c;

      CC = mdivide_left_tri_low(L_uu, K_us);
      D = mdivide_left_tri_low(B_chol, CC);
      cov = K_ss - CC'*CC + D'*D + diag_matrix(rep_vector(1e-4, n_pred));
      f_pred = multi_normal_rng(mu, cov);
    }
    return f_pred;
  }  
  
}

data {
  int<lower=1> nP;
  int<lower=1> nI;
  int<lower=1> nP_pred;
  int<lower=1> nI_pred;
  vector[nP] tP;
  vector[nI] tI;
  vector[nP_pred] tP_pred;
  vector[nI_pred] tI_pred;
  vector[nP] yP;
  vector[nI] yI;
  
  // posterior samples
  int N_samples;
  vector[N_samples] rho;
  vector[N_samples] alpha;
  real sigmaP;
  real sigmaI;
  vector[N_samples] rho_d;
  vector[N_samples] alpha_d;
  vector[N_samples] R;
  vector[N_samples] C;
  
  int<lower=1> mP;
  int<lower=1> mI;
  vector<lower=0,upper=1>[mP] zP;
  vector<lower=0,upper=1>[mI] zI;
}

parameters {
}

model {

}

generated quantities {
  vector[nP_pred] f_P;
  matrix[N_samples, nP_pred] y_P;
  vector[nI_pred] f_I;
  matrix[N_samples, nI_pred] y_I;
  
  for(n in 1:N_samples) {
    f_P = Pgp_pred_rng(tP, tI, yP, yI, tP_pred[1:nP_pred], alpha[n], rho[n], rho_d[n], alpha_d[n], R[n], C[n], sigmaP, sigmaI, zP, zI);
    for(np in 1:nP_pred) {
      y_P[n, np] = normal_rng(f_P[np], sigmaP);
    }  
  }
  
  for(n in 1:N_samples) {
    f_I = Igp_pred_rng(tP, tI, yP, yI, tI_pred[1:nI_pred], alpha[n], rho[n], rho_d[n], alpha_d[n], R[n], C[n], sigmaP, sigmaI, zP, zI);
    for(np in 1:nI_pred) {
      y_I[n, np] = normal_rng(f_I[np], sigmaI);
    }  
  }
}