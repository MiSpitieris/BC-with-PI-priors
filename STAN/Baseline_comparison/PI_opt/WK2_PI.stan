functions {
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
  // Physics informed prior kernel of the WK2 model
  matrix K_wk2(vector tP, 
               vector tI,
               real rho,
               real alpha,
               real sigmaP,
               real sigmaI,
               real R,
               real C) {
    int nP = rows(tP);
    int nI = rows(tI);
    matrix[nP + nI, nP + nI] K;

// KP
for (i in 1:(nP-1)){
  K[i,i] = pow(alpha, 0.2e1);
  for (j in (i+1):nP){
    K[i,j] = exp(-pow(tP[i] - tP[j], 0.2e1) * pow(rho, -0.2e1));
    K[i,j] = pow(alpha, 0.2e1) * K[i,j];
    K[j,i] = K[i,j];
  }
  K[nP,nP] = pow(alpha, 0.2e1);
}
K[1:nP, 1:nP] = K[1:nP, 1:nP] + diag_matrix(rep_vector(pow(sigmaP, 0.2e1), nP));

// KPI
for (i in 1:nP){
  for (j in 1:nI){
    K[i, nP + j] = 0.1e1 / R * exp(-pow(tP[i] - tI[j], 0.2e1) * pow(rho, -0.2e1)) 
    + 0.2e1 * C * (tP[i] - tI[j])* pow(rho, -0.2e1) 
    * exp(-pow(tP[i] - tI[j], 0.2e1) * pow(rho, -0.2e1));
    K[i, nP + j] = pow(alpha, 0.2e1) * K[i, nP + j];              
  }
}

         // KIP (KIP = KPI')
        // K[(nP + 1):(nP + nI), 1:nP] = K[1:nP, (nP + 1):(nP + nI)]';
        for (i in 1:nI){
          for (j in 1:nP){
            K[nP + i, j] = 0.1e1 / R * exp(-pow(tI[i] - tP[j], 0.2e1) * pow(rho, -0.2e1))
            - 0.2e1 * C * (tI[i] - tP[j]) * pow(rho, -0.2e1)
            * exp(-pow(tI[i] - tP[j], 0.2e1) * pow(rho, -0.2e1));
            K[nP + i, j] = pow(alpha, 0.2e1) * K[nP + i, j];
          }
        }        
        // KI
        for (i in 1:(nI-1)){
          K[nP + i, nP +i] = 
            pow(R, -0.2e1) * exp(-pow(tI[i] - tI[i], 0.2e1) * pow(rho, -0.2e1)) 
          + C * C * (0.2e1 * pow(rho, -0.2e1) * exp(-pow(tI[i] - tI[i], 0.2e1) 
                                                        * pow(rho, -0.2e1)) - 0.4e1 * pow(tI[i] - tI[i], 0.2e1) 
                         * pow(rho, -0.4e1) * exp(-pow(tI[i] - tI[i], 0.2e1) * pow(rho, -0.2e1)));
          K[nP + i, nP +i] = pow(alpha, 0.2e1) * K[nP + i, nP +i];
          for (j in (i+1):nI){
            K[nP + i, nP +j] = 
              pow(R, -0.2e1) * exp(-pow(tI[i] - tI[j], 0.2e1) * pow(rho, -0.2e1)) 
            + C * C * (0.2e1 * pow(rho, -0.2e1) * exp(-pow(tI[i] - tI[j], 0.2e1) 
                                                          * pow(rho, -0.2e1)) - 0.4e1 * pow(tI[i] - tI[j], 0.2e1) 
                           * pow(rho, -0.4e1) * exp(-pow(tI[i] - tI[j], 0.2e1) * pow(rho, -0.2e1)));
            K[nP + i, nP +j] = pow(alpha, 0.2e1) * K[nP + i, nP +j];
            K[nP + j, nP +i] = K[nP + i, nP +j];
          }
          K[nP + nI, nP +nI] = 
            pow(R, -0.2e1) * exp(-pow(tI[nI] - tI[nI], 0.2e1) * pow(rho, -0.2e1)) 
          + C * C * (0.2e1 * pow(rho, -0.2e1) * exp(-pow(tI[nI] - tI[nI], 0.2e1) 
                                                        * pow(rho, -0.2e1)) - 0.4e1 * pow(tI[nI] - tI[nI], 0.2e1) 
                         * pow(rho, -0.4e1) * exp(-pow(tI[nI] - tI[nI], 0.2e1) * pow(rho, -0.2e1)));
          K[nP + nI, nP +nI] = pow(alpha, 0.2e1) * K[nP + nI, nP +nI];
        }
        K[(nP + 1):(nP + nI), (nP + 1):(nP + nI)] = K[(nP + 1):(nP + nI), (nP + 1):(nP + nI)]
        + diag_matrix(rep_vector(pow(sigmaI, 0.2e1), nI));
        return cholesky_decompose(K);
  }
}

data {
  int<lower=1> nP;
  int<lower=1> nI;
  vector[nP] tP;
  vector[nI] tI;
  vector[nP] yP;
  vector[nI] yI;
}

transformed data {
  vector[nP + nI] y = append_row(yP, yI);
}

parameters {
  // hyper-parameters
  real<lower=0,upper=1> rho;
  real<lower=0> alpha;
  real<lower=0,upper=200> mu_wk2;
  real<lower=0,upper=20> sigmaP;
  real<lower=0,upper=20> sigmaI;
  // physical parameters
  real<lower=0.5, upper=3> R;
  real<lower=0.5, upper=3> C;
}

model {
  // Chol. of PI kernel
  matrix[nP + nI, nP + nI] L_K = K_wk2(tP, tI, rho, alpha, sigmaP, sigmaI, R, C);
  // mean vector
  vector[nP + nI] mu = mu_fn(mu_wk2, R, nP, nI);
  // priors
  // rho ~ normal(0,1.0/3);
  // alpha ~ normal(0,20);
  // mu_wk2 ~ normal(mean(yP), 20);
  
  y ~ multi_normal_cholesky(mu, L_K);
}
