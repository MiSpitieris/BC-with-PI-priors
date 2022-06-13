functions {
  vector mu_fn(real mu_wk2, 
               real R, 
               int nP, 
               int nI){
                 real R01 = R*2.5 + 0.5;
                 vector[nP]  mP = rep_vector(mu_wk2,nP); 
                 vector[nI]  mI = rep_vector((1/R01)*mu_wk2,nI);
                 vector[nP + nI] mu;
                 mu = append_row(mP, mI);
                 return(mu);
  }
  // Physics informed prior kernel of the WK2 model
  matrix K_wk2(matrix tP, 
               matrix tI,
               real l,
               real sigma,
               real p,
               real sigmaP,
               real sigmaQ,
               real R,
               real C) {
    int nP = rows(tP);
    int nI = rows(tI);
    matrix[nP + nI, nP + nI] K;
    real a;
    real b;
    real R01 = R*2.5 + 0.5;
    real C01 = C*2.5 + 0.5;
    // KP
    for (i in 1:(nP-1)){
      K[i,i] = pow(sigma, 0.2e1);
      for (j in (i+1):nP){
        a=tP[i,1]; 
        b=tP[j,1];
        K[i,j] = exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1));
        K[i,j] = pow(sigma, 0.2e1) * K[i,j];
        K[j,i] = K[i,j];
      }
      K[nP,nP] = pow(sigma, 0.2e1);
    }
    K[1:nP, 1:nP] = K[1:nP, 1:nP] + diag_matrix(rep_vector(pow(sigmaP, 0.2e1), nP));
    
    // KPI
    for (i in 1:nP){
      for (j in 1:nI){
        a=tP[i,1];  
        b=tI[j,1]; 
        K[i, nP + j] = 
        0.1e1 / R01 * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1)) + 
        0.4e1 * C01 * sin(pi() * (a - b) * p) * pow(l, -0.2e1) * pi() * p * 
        cos(pi() * (a - b) * p) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1));
        K[i, nP + j] = pow(sigma, 0.2e1) * K[i, nP + j];              
      }
    }
    
    // KIP (KIP = KPI')
    //K[(nP + 1):(nP + nI), 1:nP] = K[1:nP, (nP + 1):(nP + nI)]';
    for (i in 1:nI){
      for (j in 1:nP){
       a=tI[i,1];  
       b=tP[j,1]; 
       K[nP + i, j] = 
       0.1e1 / R01 * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1)) - 
       0.4e1 * C01 * sin(pi() * (a - b) * p) * pow(l, -0.2e1) * pi() * p * 
       cos(pi() * (a - b) * p) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1));
       K[nP + i, j] = pow(sigma, 0.2e1) * K[nP + i, j];
     }
   }        
   // KI
    for (i in 1:(nI-1)){
     K[nP + i, nP +i] = 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(l, -0.2e1)) * pow(R01, -0.2e1) + 
     C01 * C01 * (0.4e1 * pi() * pi() * p * p * pow(cos(pi() * 0 * p), 0.2e1) * 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(l, -0.2e1)) * pow(l, -0.2e1) - 
     0.4e1 * pow(sin(pi() * 0 * p), 0.2e1) * pi() * pi() * p * p * 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(l, -0.2e1)) * pow(l, -0.2e1) - 
     0.16e2 * pow(sin(pi() * 0 * p), 0.2e1) * pi() * pi() * p * p * 
     pow(cos(pi() * 0 * p), 0.2e1) * exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * 
     pow(l, -0.2e1)) * pow(l, -0.4e1));
      K[nP + i, nP +i] = pow(sigma, 0.2e1) * K[nP + i, nP +i];
     // if(yI[i]!=yDia) K[nP + i, nP +i] = K[nP + i, nP +i] + pow(sigmaQ,0.2e1);
     // if(yI[i]==yDia) K[nP + i, nP +i] = K[nP + i, nP +i] + 1e-12;          
     for (j in (i+1):nI){
      a=tI[i,1];
      b=tI[j,1];
      K[nP + i, nP +j] = 
      exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1)) * pow(R01, -0.2e1) + 
      C01 * C01 * (0.4e1 * pi() * pi() * p * p * pow(cos(pi() * (a - b) * p), 0.2e1) * 
      exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1)) * pow(l, -0.2e1) - 
      0.4e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pi() * pi() * p * p * 
      exp(-0.2e1 * pow(sin(pi() * (a - b) *       p), 0.2e1) * pow(l, -0.2e1)) * 
      pow(l, -0.2e1) - 0.16e2 * pow(sin(pi() * (a - b) * p), 0.2e1) * pi() * pi() * p * p * 
      pow(cos(pi() * (a - b) * p), 0.2e1) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * 
      pow(l, -0.2e1)) * pow(l, -0.4e1));

      K[nP + i, nP +j] = pow(sigma, 0.2e1) * K[nP + i, nP +j];
      K[nP + j, nP +i] = K[nP + i, nP +j];
     }
     K[nP + nI, nP +nI] = 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(l, -0.2e1)) * pow(R01, -0.2e1) + 
     C01 * C01 * (0.4e1 * pi() * pi() * p * p * pow(cos(pi() * 0 * p), 0.2e1) * 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(l, -0.2e1)) * pow(l, -0.2e1) - 
     0.4e1 * pow(sin(pi() * 0 * p), 0.2e1) * pi() * pi() * p * p * 
     exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * pow(l, -0.2e1)) * pow(l, -0.2e1) - 
     0.16e2 * pow(sin(pi() * 0 * p), 0.2e1) * pi() * pi() * p * p * 
     pow(cos(pi() * 0 * p), 0.2e1) * exp(-0.2e1 * pow(sin(pi() * 0 * p), 0.2e1) * 
     pow(l, -0.2e1)) * pow(l, -0.4e1));
     K[nP + nI, nP +nI] = pow(sigma, 0.2e1) * K[nP + nI, nP +nI];
     // if(yI[nI]!=yDia) K[nP + nI, nP +nI] = K[nP + nI, nP +nI] + pow(sigmaQ,0.2e1);
     // if(yI[nI]==yDia) K[nP + nI, nP +nI] = K[nP + nI, nP +nI] + 1e-12;
     }
     K[(nP + 1):(nP + nI), (nP + 1):(nP + nI)] = K[(nP + 1):(nP + nI), (nP + 1):(nP + nI)]
            + diag_matrix(rep_vector(pow(sigmaQ, 0.2e1), nI));
     return cholesky_decompose(K);
  }
  
  matrix KP_pred(matrix tP,
                 matrix tP_pred,
                 real l,
                 real sigma,
                 real p){
    int nP = rows(tP);
    int nP_pred = rows(tP_pred);
    real a;
    real b;
    matrix[nP_pred, nP] KP;
    
    for (i in 1:nP_pred){
      for (j in 1:nP){
        a=tP_pred[i,1];
        b=tP[j,1];
        KP[i,j] = exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1));
        KP[i,j] = pow(sigma, 0.2e1) * KP[i,j];
      }
    }
    return KP;
  }
  
  // KPI_pred
  matrix KPI_pred(matrix tI,
                  matrix tP_pred,
                  real l,
                  real sigma,
                  real p,
                  real R,
                  real C
  ){
    int nP = rows(tP_pred);
    matrix[nP,1] tP = tP_pred;
    int nI = rows(tI);
    matrix[nP, nI] KPI;
    real a;
    real b;
    real R01 = R*2.5 + 0.5;
    real C01 = C*2.5 + 0.5;
    
    for (i in 1:nP){
      for (j in 1:nI){
        a=tP[i,1]; 
        b=tI[j,1];
        KPI[i, j] = 
        0.1e1 / R01 * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1)) + 
        0.4e1 * C01 * sin(pi() * (a - b) * p) * pow(l, -0.2e1) * pi() * p * 
        cos(pi() * (a - b) * p) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1));
        KPI[i, j] = pow(sigma, 0.2e1) * KPI[i, j];
      }
    }
    return KPI; 
  }
  
  // KIP_pred
  matrix KIP_pred(matrix tI_pred,
                  matrix tP,
                  real l,
                  real sigma,
                  real p,
                  real R,
                  real C){
    int nI = rows(tI_pred);
    matrix[nI,1] tI = tI_pred;
    int nP = rows(tP);
    matrix[nI, nP] KIP;
    real a;
    real b;
    real R01 = R*2.5 + 0.5;
    real C01 = C*2.5 + 0.5;
    for (i in 1:nI){
      for (j in 1:nP){
        a=tI[i,1]; 
        b=tP[j,1];
        KIP[i, j] = 
        0.1e1 / R01 * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1)) - 
        0.4e1 * C01 * sin(pi() * (a - b) * p) * pow(l, -0.2e1) * pi() * p * 
        cos(pi() * (a - b) * p) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1));
        KIP[i, j] = pow(sigma, 0.2e1) * KIP[i, j];
      }
    }
    return KIP;
  }
  
  matrix KI_pred(matrix tI_pred,
                 matrix tI,
                 real l,
                 real sigma,
                 real p,
                 real R,
                 real C){
    int nI = rows(tI);
    int nI_pred = rows(tI_pred);
    matrix[nI_pred, nI] KI;
    real a;
    real b;
    real R01 = R*2.5 + 0.5;
    real C01 = C*2.5 + 0.5;
    // KI
    for (i in 1:nI_pred){
      for (j in 1:nI){
        a=tI_pred[i,1]; 
        b=tI[j,1];
        KI[i, j] = 
        exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1)) * pow(R01, -0.2e1) + 
        C01 * C01 * (0.4e1 * pi() * pi() * p * p * pow(cos(pi() * (a - b) * p), 0.2e1) * 
        exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pow(l, -0.2e1)) * pow(l, -0.2e1) - 
        0.4e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * pi() * pi() * p * p * 
        exp(-0.2e1 * pow(sin(pi() * (a - b) *       p), 0.2e1) * pow(l, -0.2e1)) * 
        pow(l, -0.2e1) - 0.16e2 * pow(sin(pi() * (a - b) * p), 0.2e1) * pi() * pi() * p * p * 
        pow(cos(pi() * (a - b) * p), 0.2e1) * exp(-0.2e1 * pow(sin(pi() * (a - b) * p), 0.2e1) * 
        pow(l, -0.2e1)) * pow(l, -0.4e1));
        KI[i, j] = pow(sigma, 0.2e1) * KI[i, j];
      }
    }
    return KI;
  }
  
  vector Pgp_pred_rng(matrix tP,
                      matrix tI,
                      vector yP,
                      vector yI,
                      matrix tP_pred,
                      real l,
                      real sigma,
                      real p,
                      real R,
                      real C,
                      real sigmaP,
                      real sigmaQ) {
    int nP = rows(tP);
    int nI = rows(tI);
    int N = nP + nI;
    int nP_pred = rows(tP_pred);
    vector[nP_pred] f2;
    {
      matrix[N, N] L_K;
      vector[N] K_div_y;
      matrix[nP_pred, N] qT_p;
      matrix[N, nP_pred] v_pred;
      vector[nP_pred] f2_mu;
      matrix[nP_pred, nP_pred] cov_f2;
      matrix[nP_pred, nP_pred] diag_delta;
      vector[N] y;
      y[1:nP] = yP[1:nP];
      y[(nP+1):N] = yI[1:nI];
       
      L_K =  K_wk2(tP, tI, l, sigma, p, sigmaP, sigmaQ, R, C);
      K_div_y = mdivide_left_tri_low(L_K, y);
      K_div_y = mdivide_right_tri_low(K_div_y', L_K)';

      qT_p[1:nP_pred, 1:nP] = KP_pred(tP, tP_pred, l, sigma, p);
      qT_p[1:nP_pred, (nP+1):N] = KPI_pred(tI, tP_pred, l, sigma, p, R, C);
      f2_mu = (qT_p * K_div_y);
      v_pred = mdivide_left_tri_low(L_K, qT_p');
      cov_f2 = KP_pred(tP_pred, tP_pred, l, sigma, p) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(1e-6, nP_pred));
                                                                      
      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
  
  vector Igp_pred_rng(matrix tP,
                      matrix tI,
                      vector yP,
                      vector yI,
                      matrix tI_pred,
                      real l,
                      real sigma,
                      real p,
                      real R,
                      real C,
                      real sigmaP,
                      real sigmaQ) {
    int nP = rows(tP);
    int nI = rows(tI);
    int N = nP + nI;
    int nI_pred = rows(tI_pred);
    vector[nI_pred] f2;
    {
      matrix[N, N] L_K;
      vector[N] K_div_y;
      matrix[nI_pred, N] qT_I;
      matrix[N, nI_pred] v_pred;
      vector[nI_pred] f2_mu;
      matrix[nI_pred, nI_pred] cov_f2;
      matrix[nI_pred, nI_pred] diag_delta;
      vector[N] y;

      y[1:nP] = yP[1:nP];
      y[(nP+1):N] = yI[1:nI];
      
      L_K =  K_wk2(tP, tI, l, sigma, p, sigmaP, sigmaQ, R, C);
      K_div_y = mdivide_left_tri_low(L_K, y);
      K_div_y = mdivide_right_tri_low(K_div_y', L_K)';
      qT_I[1:nI_pred, 1:nP] = KIP_pred(tI_pred, tP, l, sigma, p, R, C);
      qT_I[1:nI_pred, (nP+1):N] = KI_pred(tI_pred, tI, l, sigma, p, R, C);
      f2_mu = (qT_I * K_div_y);
      v_pred = mdivide_left_tri_low(L_K, qT_I');
      cov_f2 = KI_pred(tI_pred, tI_pred, l, sigma, p, R, C) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(1e-6, nI_pred));
                                                                    
      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
}

data {
  int<lower=1> nP;
  int<lower=1> nI;
  int<lower=1> nP_pred;
  int<lower=1> nI_pred;
  matrix[nP,1] tP;
  matrix[nI,1] tI;
  matrix[nP_pred,1] tP_pred;
  matrix[nI_pred,1] tI_pred;
  vector[nP] yP;
  vector[nI] yI;
}

transformed data {
  vector[nP + nI] y = append_row(yP, yI);
}

parameters {
  real<lower=0> l;
  real<lower=0,upper=1> sigma;
  real<lower=0.8, upper=1.2> p;
  real<lower=0> sigmaP;
  real<lower=0> sigmaQ;
  real<lower=0> mu_wk2;
  real<lower=0, upper=1> R;
  real<lower=0, upper=1> C;
}

model {
  // Chol. of PI kernel
  matrix[nP + nI, nP + nI] L_K = K_wk2(tP, tI, l, sigma, p, sigmaP, sigmaQ, R, C);
  // mean vector
  vector[nP + nI] mu = mu_fn(mu_wk2, R, nP, nI);
  // priors
  l~normal(0,1);
  sigma ~ normal(0,1.0/3);
  sigmaP ~ normal(0, 0.1/3);
  sigmaQ ~ normal(0, 0.1/3);
  mu_wk2 ~ normal(mean(yP), 0.3);
  
  y ~ multi_normal_cholesky(mu, L_K);
}

generated quantities {
  vector[nP_pred] f_P;
  vector[nP_pred] y_P;
  vector[nI_pred] f_I;
  vector[nI_pred] y_I;

  f_P = Pgp_pred_rng(tP, tI, yP, yI, tP_pred, l, sigma, p, R, C, sigmaP, sigmaQ);
  for (n1 in 1:nP_pred)
    y_P[n1] = normal_rng(f_P[n1], sigmaP);

  f_I = Igp_pred_rng(tP, tI, yP, yI, tI_pred, l, sigma, p, R, C, sigmaP, sigmaQ);
  for (n2 in 1:nI_pred)
    y_I[n2] = normal_rng(f_I[n2], sigmaQ);
}

