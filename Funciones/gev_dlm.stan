data {
  int<lower=0> n;
  vector[n] y;
  real m0;
  real phi;
  real<lower=0> C0;
}
parameters {
  real k;
  real nu0;
  real<lower=0> tau;
  real<lower=0> scale;
  real<lower=0> sigma;
  vector[n] eta;
  vector[n] nu;
}
transformed parameters{
  vector[n] g_t;
  vector[n+1] mu;
  
  mu[1] = m0 + C0*nu0;
  
  for(i in 1:n){
     mu[i+1] = phi*mu[i] + tau*nu[i];
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
  target += normal_lpdf(nu | 0, 1);
  
  // likelihood
  target += normal_lpdf(y | g_t, sigma);
}
