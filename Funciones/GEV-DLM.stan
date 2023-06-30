data {
  int<lower=0> n;
  int<lower=0> d;
  vector[n] y;
  vector[d] FF;
  vector[d] m0;
  matrix[d,d] G;
  matrix[d,d] C0;
}
parameters {
  real k;
  real<lower=0> scale;
  real<lower=0> sigma;
  
  vector [d] nu0;
  vector [n] eta;
  vector<lower=0> [d] tau;
  array[n] vector [d] nu;
}
transformed parameters{
  vector[n] g_t;
  vector[n] mu;
  array[n+1] vector [d] mt;
  matrix[d,d] W = diag_matrix(pow(tau,2));
  
  mt[1] = m0 + C0*nu0;
  
  for(i in 1:n){
     mt[i+1] = G*mt[i] + W*nu[i];
     
     mu[i] = sum(FF.* mt[i]);
     g_t[i] = mu[i] + scale*pow(k,-1)*(exp(k*eta[i]) - 1);
  }
}
model {
  // priors
  target += normal_lpdf(k | 0, 1);
  target += normal_lpdf(nu0 | 0, 1);
  // Variances
  target += student_t_lpdf(tau | 3, 0, 1);
  target += student_t_lpdf(scale | 3, 0, 1);
  target += normal_lpdf(sigma | 0, 0.5);
  
  // latent varaibles
  target += gumbel_lpdf(eta | 0, 1);
  for(i in 1:n) target += normal_lpdf(nu[i] | 0, 1);
  
  // likelihood
  target += normal_lpdf(y | g_t, sigma);
}

