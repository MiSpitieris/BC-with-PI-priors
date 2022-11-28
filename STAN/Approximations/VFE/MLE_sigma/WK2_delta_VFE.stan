// Variational Free Energy (VFE) GP || Titsias 2009
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
     + C2 * (0.2e1 * rho2_ * exp(-pow(tI[i] - tI[i], 0.2e1)
                                                  * rho2_) - 0.4e1 * pow(tI[i] - tI[i], 0.2e1)
                    * rho4_ * exp(-pow(tI[i] - tI[i], 0.2e1) * rho2_));
     K[nP + i, nP +i] = alpha2 * K[nP + i, nP +i];
     for (j in (i+1):nI){
       K[nP + i, nP +j] =
       R2_ * exp(-pow(tI[i] - tI[j], 0.2e1) * rho2_)
       + C2 * (0.2e1 * rho2_ * exp(-pow(tI[i] - tI[j], 0.2e1)
                                                          * rho2_) - 0.4e1 * pow(tI[i] - tI[j], 0.2e1)
                           * rho4_ * exp(-pow(tI[i] - tI[j], 0.2e1) * rho2_));
     K[nP + i, nP +j] = alpha2 * K[nP + i, nP +j];
     K[nP + j, nP +i] = K[nP + i, nP +j];
     }
     K[nP + nI, nP +nI] =
      R2_ * exp(-pow(tI[nI] - tI[nI], 0.2e1) * rho2_)
      + C2 * (0.2e1 * rho2_ * exp(-pow(tI[nI] - tI[nI], 0.2e1)
                                                        * rho2_) - 0.4e1 * pow(tI[nI] - tI[nI], 0.2e1)
                         * rho4_ * exp(-pow(tI[nI] - tI[nI], 0.2e1) * rho2_));
      K[nP + nI, nP +nI] = alpha2 * K[nP + nI, nP +nI];
    }
    return K;
  }
  
  matrix K_diag(vector tP,
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
    matrix[nP + nI, nP + nI] K;

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
  K = diag_matrix(v);
  return K;
 }
}
data {
  int<lower=1> nP;
  int<lower=1> nI;
  int<lower=1> mP;
  int<lower=1> mI;
  vector[nP] tP;
  vector[nI] tI;
  vector[nP] yP;
  vector[nI] yI;
  
  real sigmaP;
  real sigmaI;
}
transformed data {
  vector[nP + nI] y = append_row(yP, yI);
  int N = nP+nI;
  int m = mP+mI;
}
parameters {
  // hyper-parameters
  real<lower=0.05> rho;
  real<lower=0,upper=100> alpha;
  real<lower=0.05,upper=1> rho_d;
  real<lower=0,upper=50> alpha_d;
  real<lower=0,upper=400> mu_wk2;
  // physical parameters
  real<lower=0.5, upper=3> R;
  real<lower=0.5, upper=3> C;
  // ind
  vector<lower=0,upper=1>[mP] zP;
  vector<lower=0,upper=1>[mI] zI;
}

model {
  matrix[N,m] K_fu;
  matrix[m,N] K_uf;
  matrix[N,N] K_ff;
  matrix[m,m] K_uu;
  matrix[m,m] L_uu; 
  matrix[N,N] Sigma_inv; 
  matrix[m,N] A;
  matrix[m,m] B;
  matrix[m,m] B_chol;
  matrix[m,m] AS_invA;
  vector[N] sigma_vec;
  vector[N] sigma_vec_inv;
  real log_pi;
  real log_det_B;
  real log_det_noise;
  vector[N] r;
  vector[N] r_l;
  vector[m] c;
  real quad;
  real trc;
  
  K_fu = K_wk2_XZ(tP, tI, zP, zI, rho, alpha, rho_d, alpha_d, R, C);
  K_uf = K_fu';
  K_uu = K_wk2(zP, zI, rho, alpha, rho_d, alpha_d, R, C)+ diag_matrix(rep_vector(1e-6,m));
  K_ff = K_diag(tP, tI, rho, alpha, rho_d, alpha_d, R, C);
  
  sigma_vec = to_vector(append_row(rep_vector(square(sigmaP),nP),rep_vector(square(sigmaI),nI)));
  sigma_vec_inv = to_vector(append_row(rep_vector(1/square(sigmaP),nP),rep_vector(1/square(sigmaI),nI)));
  Sigma_inv = diag_matrix(sigma_vec_inv);
  
  r = y- mu_fn(mu_wk2, R, nP, nI);
  r_l = Sigma_inv*r;
  
  L_uu = cholesky_decompose(K_uu);
  A = mdivide_left_tri_low(L_uu, K_uf);
  AS_invA = A*Sigma_inv*A';
  B = diag_matrix(rep_vector(1,m)) + AS_invA;
  B_chol = cholesky_decompose(B+ diag_matrix(rep_vector(1e-6,m)));
  c = mdivide_left_tri_low(B_chol, A*r_l);
  
  quad = -0.5*(r'*r_l -c'*c);
  
  // det terms
  log_pi = -(N/2.0)*log(2*pi());
  log_det_B = - 0.5*sum(log(diagonal(B_chol)));
  log_det_noise = -0.5*(nP * log(square(sigmaP))+ nI * log(square(sigmaI)));
  
  // trace term
  trc = -0.5*(trace(Sigma_inv*K_ff)+ trace(AS_invA)); 
   
  // priors
  rho ~ normal(0,1.0/3);
  alpha ~ normal(0,20);
  rho_d ~ normal(0,1.0/3);
  alpha_d ~ normal(0,20);
  mu_wk2 ~ normal(mean(yP), 20);

  target += log_pi + log_det_B +log_det_noise + quad +trc;
}

