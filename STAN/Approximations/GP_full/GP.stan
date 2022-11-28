// Sparse GP Stan Code
data {
  int<lower=1> N;
  real x[N];
  vector[N] y;
}

parameters {
  real<lower=0,upper=1> rho;
  real<lower=0, upper=100> alpha;
  real<lower=0,upper=20> sigma;
}

model {
  matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)+ diag_matrix(rep_vector(square(sigma), N));
  matrix[N, N] L_cov = cholesky_decompose(cov);

  y ~ multi_normal_cholesky(rep_vector(0, N), L_cov);
}