data {
  int<lower=0> n;
  vector[n] y;
}
parameters {
  real k;
  real mu;
  real<lower=0> scale;
  real<lower=0> sigma;
  vector[n] eta;
}
model {
  vector[n] g_t;
   for(i in 1:n){
     g_t[i] = mu + scale*pow(k,-1)*(exp(k*eta[i]) - 1);
  }
  // priors
  target += normal_lpdf(k | 0, 1);
  target += normal_lpdf(mu | 0, 1);
  // Variances
  target += student_t_lpdf(scale | 3, 0, 1);
  target += normal_lpdf(sigma | 0, 0.5);
  
  // latent varaibles
  target += gumbel_lpdf(eta | 0, 1);
  
  // likelihood
  target += normal_lpdf(y | g_t, sigma);
}
