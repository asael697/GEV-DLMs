
data {
  int<lower=0> n;
  int<lower=0> d;
  int<lower=0> k;
  array[n] vector [d] y;
  matrix[k,d] FF;
  matrix[k,k] G;
  vector[k] m0;
  matrix[k,k] C0;
}
parameters {
  real<lower=0> tau;
  real<lower=0> sigma;
}
transformed parameters{
  matrix[k,k] W; 
  matrix[d,d] V; 
  
  // DLM Transcitions
  vector[k] a;
  vector[d] error;
    
  //array of vectors for the mean
  array[n+1] vector [k] mt;
  array[n] vector [d] ft;
  
  matrix[k,d] A;
  matrix[k,k] R;
    
  array[n+1] matrix[k,k] Ct;
  array[n]   matrix[d,d] Qt;
  
  W =  pow(tau,2)*diag_matrix(rep_vector(1,k));
  V = pow(sigma,2)*diag_matrix(rep_vector(1,d));
  
  mt[1] = m0;
  Ct[1] = C0;
  
  for(i in 1:n){
    // prior at t
    a = G * mt[i];
    R = (G *Ct[i] *G') + W;
    
    // Marginal likelihood: 
    ft[i] = FF' * a;
    Qt[i] = (FF' * R * FF) + V;
    
    // additional computations at t
    error = y[i] - ft[i];
    A = R * FF / Qt[i];
    
    //  Posterior at t
    mt[i+1] = a + (A*error);
    Ct[i+1] = R - (A  * Qt[i] * A');
  }
}
model {
  // priors
  target += normal_lpdf(sigma | 0, 1);
  target += normal_lpdf(tau   | 0, 1);
  
  // likelihood
   for (i in 1:n) {
     target += multi_normal_lpdf(y[i] | ft[i], Qt[i]);
   }
}
generated quantities{
  array[n+1] vector[k] mu;
  array[n] vector[d] y_pred;
  matrix[k,k] At;
  matrix[k,k] Q;
  matrix[k,k] H;
  vector[k] h;

  mu[n+1] = multi_normal_rng(mt[n+1],Ct[n+1]);
 
  for(i in 1:n){
    At = Ct[n-i+1] * G';
    Q = inverse(G*At + W);
    H = Ct[n-i+1] - At*Q*At';
    
    h = mt[n-i+1] + At*Q*(mu[n-i+2] - G*mt[n-i+1]);
    mu[n-i+1] = multi_normal_rng(h,H);
    y_pred[i] = multi_normal_rng(ft[i],Qt[i]);
  }
}
