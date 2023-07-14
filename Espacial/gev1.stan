functions{
  real gev_lpdf(vector y, real mu, real sigma, real xi) {
    int n = num_elements(y);
    vector[n] z;
    vector[n] logf;
    
    for(i in 1:n){
      z[i] = pow(sigma,-1)*(y[i] -mu);
      z[i] = xi==0 ? exp(z[i]) : pow(1 + (xi *z[i]), -1/xi);
      logf[i] = -log(sigma) + (xi + 1) * log(z[i]) - z[i];
    }
    return sum(logf);
  }
}
data {
  int<lower=0> n;
  vector[n] y;
}
parameters {
  real k;
  real<lower=0> scale;
}
model {
  // priors
  target += normal_lpdf(k | 0, 0.5);
  // Variances
  target += student_t_lpdf(scale | 3, 0, 1);

  // Data Model
  target += gev_lpdf(y | 0, scale, k);
}
