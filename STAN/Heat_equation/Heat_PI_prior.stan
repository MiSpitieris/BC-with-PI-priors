functions {
  // mean function
  vector mu_fn(real mu_heat, 
               int nU, 
               int nF){
                 vector[nU]  mU = rep_vector(mu_heat,nU); 
                 vector[nF]  mF = rep_vector(0,nF);
                 vector[nU + nF] mn;
                 mn = append_row(mU, mF);
                 return(mn);
  }
  // K11
  real Kuu_fn(real x1, real x2,
             real t1, real t2,
             real lx, real lt){
               real val;
               real dx = x1-x2;
               real dt = t1-t2;
               val = exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1);
               return(val);
  }
  // K12
  real Kuf_fn(real x1, real x2,
              real t1, real t2,
              real lx, real lt,
              real k){
                real val;
                real dx = x1-x2;
                real dt = t1-t2;
                val = exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * (dt) * pow(lt, -0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1) - 
                      k * (-pow(lx, -0.2e1) * exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1) + 
                      pow(dx, 0.2e1) * pow(lx, -0.4e1) * exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1));
                return(val);
  }  
  // K21
  real Kfu_fn(real x1, real x2,
              real t1, real t2,
              real lx, real lt,
              real k){
                real val;
                real dx = x1-x2;
                real dt = t1-t2;
                val = -exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * (dt) * pow(lt, -0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1) - 
                      k * (-pow(lx, -0.2e1) * exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1) + pow(dx, 0.2e1) * 
                      pow(lx, -0.4e1) * exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1));
                return(val);
  }   
  // K22
  real Kff_fn(real x1, real x2,
              real t1, real t2,
              real lx, real lt,
              real k){
                real val;
                real dx = x1-x2;
                real dt = t1-t2;
                val = exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * pow(lt, -0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1) - 
                      exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * pow(dt, 0.2e1) * pow(lt, -0.4e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1) + 
                      k * k * (0.3e1 * pow(lx, -0.4e1) * exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1) - 
                      0.6e1 * pow(lx, -0.6e1) * pow(dx, 0.2e1) * exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1) + 
                      pow(dx, 0.4e1) * pow(lx, -0.8e1) * exp(-pow(dx, 0.2e1) * pow(lx, -0.2e1) / 0.2e1) * exp(-pow(dt, 0.2e1) * pow(lt, -0.2e1) / 0.2e1));
                return(val);
  }   
             
  // Physics informed prior kernel of the Heat equation
  matrix K_heat(matrix xtU, matrix xtF,
               real lx, real lt,
               real sigma,
               real sigmaU,
               real sigmaF,
               real alpha) {
    int nU = rows(xtU);
    int nF = rows(xtF);
    matrix[nU + nF, nU + nF] K;
    
    // KU
    for (i in 1:(nU-1)){
      K[i,i] = pow(sigma, 0.2e1);
      for (j in (i+1):nU){
        K[i,j] = Kuu_fn(xtU[i,1], xtU[j,1], xtU[i,2], xtU[j,2], lx, lt);
        K[i,j] = pow(sigma, 0.2e1) * K[i,j];
        K[j,i] = K[i,j];
      }
      K[nU,nU] = pow(sigma, 0.2e1);
    }
    K[1:nU, 1:nU] = K[1:nU, 1:nU] + diag_matrix(rep_vector(pow(sigmaU, 0.2e1), nU));
    
    // KUF
    for (i in 1:nU){
      for (j in 1:nF){
        K[i, nU + j] = Kuf_fn(xtU[i,1], xtF[j,1], xtU[i,2], xtF[j,2], lx, lt, alpha);
        K[i, nU + j] = pow(sigma, 0.2e1) * K[i, nU + j];              
      }
    }
    
    // KFU (KFU = KUF')
    K[(nU + 1):(nU + nF), 1:nU] = K[1:nU, (nU + 1):(nU + nF)]';
   //  for (i in 1:nF){
   //    for (j in 1:nU){
   //      K[nU + i, j] = Kfu_fn(xtF[i,1], xtU[j,1], xtF[i,2], xtU[j,2], lx, lt, alpha);
   //      K[nU + i, j] = pow(sigma, 0.2e1) * K[nU + i, j];
   //   }
   // }
   // KF
    for (i in 1:(nF-1)){
     K[nU + i, nU +i] = Kff_fn(xtF[i,1], xtF[i,1], xtF[i,2], xtF[i,2], lx, lt, alpha);
     K[nU + i, nU +i] = pow(sigma, 0.2e1) * K[nU + i, nU +i];
     for (j in (i+1):nF){
      K[nU + i, nU +j] = Kff_fn(xtF[i,1], xtF[j,1], xtF[i,2], xtF[j,2], lx, lt, alpha);
      K[nU + i, nU +j] = pow(sigma, 0.2e1) * K[nU + i, nU +j];
      K[nU + j, nU +i] = K[nU + i, nU +j];
     }
     K[nU + nF, nU +nF] = Kff_fn(xtF[nF,1], xtF[nF,1], xtF[nF,2], xtF[nF,2], lx, lt, alpha);
     K[nU + nF, nU +nF] = pow(sigma, 0.2e1) * K[nU + nF, nU +nF];
     }
     K[(nU + 1):(nU + nF), (nU + 1):(nU + nF)] = K[(nU + 1):(nU + nF), (nU + 1):(nU + nF)]
            + diag_matrix(rep_vector(pow(sigmaF, 0.2e1), nF));
     return cholesky_decompose(K);
  }
  
  matrix KU_pred(matrix xtU,
                 matrix xtU_pred,
                 real lx, real lt,
                 real sigma
  ){
    int nU = rows(xtU);
    int nU_pred = rows(xtU_pred);
    matrix[nU_pred, nU] KU;
    
    for (i in 1:nU_pred){
      for (j in 1:nU){
        KU[i,j] = Kuu_fn(xtU_pred[i,1], xtU[j,1], xtU_pred[i,2], xtU[j,2], lx, lt);
        KU[i,j] = pow(sigma, 0.2e1) * KU[i,j];
      }
    }
    return KU;
  }
  
  matrix KUF_pred(matrix xtF,
                  matrix xtU_pred,
                  real lx, real lt,
                  real sigma,
                  real alpha
  ){
    int nU = rows(xtU_pred);
    matrix[nU,2] xtU = xtU_pred;
    int nF = rows(xtF);
    matrix[nU, nF] KUF;
    
    // KUF
    for (i in 1:nU){
      for (j in 1:nF){
        KUF[i, j] = Kuf_fn(xtU[i,1], xtF[j,1], xtU[i,2], xtF[j,2], lx, lt, alpha);
        KUF[i, j] = pow(sigma, 0.2e1) * KUF[i, j];              
      }
    }
    return KUF; 
  }
  
  matrix KFU_pred(matrix xtF_pred,
                  matrix xtU,
                  real lx, real lt,
                  real sigma,
                  real alpha){
    int nF = rows(xtF_pred);
    matrix[nF,2] xtF = xtF_pred;
    int nU = rows(xtU);
    matrix[nF, nU] KFU;

    for (i in 1:nF){
      for (j in 1:nU){
        KFU[i, j] = Kfu_fn(xtF[i,1], xtU[j,1], xtF[i,2], xtU[j,2], lx, lt, alpha);
        KFU[i, j] = pow(sigma, 0.2e1) * KFU[i, j];
      }
    }
    return KFU;
  }
  
  matrix KF_pred(matrix xtF_pred,
                 matrix xtF,
                 real lx, real lt,
                 real sigma,
                 real alpha){
    int nF = rows(xtF);
    int nF_pred = rows(xtF_pred);
    matrix[nF_pred, nF] KF;

    // KI
    for (i in 1:nF_pred){
      for (j in 1:nF){
        KF[i, j] = Kff_fn(xtF_pred[i,1], xtF[j,1], xtF_pred[i,2], xtF[j,2], lx, lt, alpha);
        KF[i, j] = pow(sigma, 0.2e1) * KF[i, j];
      }
    }
    return KF;
  }
  
  vector Ugp_pred_rng(matrix xtU,
                      matrix xtF,
                      vector yU,
                      vector yF,
                      matrix xtU_pred,
                      real sigma,
                      real lx, 
                      real lt,
                      real alpha,
                      real sigmaU,
                      real sigmaF) {
    int nU = rows(xtU);
    int nF = rows(xtF);
    int N = nU + nF;
    int nU_pred = rows(xtU_pred);
    vector[nU_pred] f2;
    {
      matrix[N, N] L_K;
      vector[N] K_div_y;
      matrix[nU_pred, N] qT_p;
      matrix[N, nU_pred] v_pred;
      vector[nU_pred] f2_mu;
      matrix[nU_pred, nU_pred] cov_f2;
      matrix[nU_pred, nU_pred] diag_delta;
      vector[N] y;
      y[1:nU] = yU[1:nU];
      y[(nU+1):N] = yF[1:nF];
       
      L_K =  K_heat(xtU, xtF, lx, lt, sigma, sigmaU, sigmaF, alpha);
      K_div_y = mdivide_left_tri_low(L_K, y);
      K_div_y = mdivide_right_tri_low(K_div_y', L_K)';

      qT_p[1:nU_pred, 1:nU] = KU_pred(xtU, xtU_pred, lx, lt, sigma);
      qT_p[1:nU_pred, (nU+1):N] = KUF_pred(xtF, xtU_pred, lx, lt, sigma, alpha);
      f2_mu = (qT_p * K_div_y);
      v_pred = mdivide_left_tri_low(L_K, qT_p');
      
      cov_f2 = KU_pred(xtU_pred, xtU_pred, lx, lt, sigma) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(1e-6, nU_pred));
                                                                      
      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
  
  vector Fgp_pred_rng(matrix xtU,
                      matrix xtF,
                      vector yU,
                      vector yF,
                      matrix xtF_pred,
                      real sigma,
                      real lx, 
                      real lt,
                      real alpha,
                      real sigmaU,
                      real sigmaF) {
    int nU = rows(xtU);
    int nF = rows(xtF);
    int N = nU + nF;
    int nF_pred = rows(xtF_pred);
    vector[nF_pred] f2;
    {
      matrix[N, N] L_K;
      vector[N] K_div_y;
      matrix[nF_pred, N] qT_I;
      matrix[N, nF_pred] v_pred;
      vector[nF_pred] f2_mu;
      matrix[nF_pred, nF_pred] cov_f2;
      matrix[nF_pred, nF_pred] diag_delta;
      vector[N] y;

      y[1:nU] = yU[1:nU];
      y[(nU+1):N] = yF[1:nF];
      
      L_K =  K_heat(xtU, xtF, lx, lt, sigma, sigmaU, sigmaF, alpha);
      K_div_y = mdivide_left_tri_low(L_K, y);
      K_div_y = mdivide_right_tri_low(K_div_y', L_K)';
      
      qT_I[1:nF_pred, 1:nU] = KFU_pred(xtF_pred, xtU, lx, lt, sigma, alpha);
      
      qT_I[1:nF_pred, (nU+1):N] = KF_pred(xtF_pred, xtF, lx, lt, sigma, alpha);
      f2_mu = (qT_I * K_div_y);
      v_pred = mdivide_left_tri_low(L_K, qT_I');
      cov_f2 = KF_pred(xtF_pred, xtF_pred, lx, lt, sigma, alpha) - v_pred' * v_pred;
      diag_delta = diag_matrix(rep_vector(1e-6, nF_pred));
                                                                    
      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_delta);
    }
    return f2;
  }
}

data {
  int<lower=1> nU;
  int<lower=1> nF;
  int<lower=1> nU_pred;
  int<lower=1> nF_pred;
  matrix[nU,2] xtU;
  matrix[nF,2] xtF;
  matrix[nU_pred,2] xtU_pred;
  matrix[nF_pred,2] xtF_pred;
  vector[nU] yU;
  vector[nF] yF;
}

transformed data {
  vector[nU + nF] y = append_row(yU, yF);
}

parameters {
  // hyper-parameters
  real<lower=0,upper=1> lx;
  real<lower=0,upper=10> lt;
  real<lower=0,upper=100> sigma;
  real<lower=-100,upper=100> mu_heat;
  real<lower=0, upper=3> sigmaU;
  real<lower=0, upper=3> sigmaF;
  // physical parameters
  real<lower=0.01, upper=10> alpha;
}

model {
  // Chol. of PI kernel
  matrix[nU + nF, nU + nF] L_K = K_heat(xtU, xtF, lx, lt, sigma, sigmaU, sigmaF, alpha);
  // mean vector
  vector[nU + nF] mu = mu_fn(mu_heat, nU, nF);
  // priors
  lx ~ normal(0,1.0/3);
  lt ~ normal(0,1.0);
  sigma ~ normal(0,1.0/3);
  mu_heat ~ normal(mean(yU), 1);
  
  y ~ multi_normal_cholesky(mu, L_K);
}

generated quantities {
  vector[nU_pred] f_U;
  vector[nU_pred] y_U;
  vector[nF_pred] f_F;
  vector[nF_pred] y_F;

  f_U = Ugp_pred_rng(xtU, xtF, yU, yF, xtU_pred, sigma, lx, lt, alpha, sigmaU, sigmaF);
  for (n1 in 1:nU_pred)
    y_U[n1] = normal_rng(f_U[n1], sigmaU);
  
  f_F = Fgp_pred_rng(xtU, xtF, yU, yF, xtF_pred, sigma, lx, lt, alpha, sigmaU, sigmaF);
  for (n2 in 1:nF_pred)
    y_F[n2] = normal_rng(f_F[n2], sigmaF);
}
