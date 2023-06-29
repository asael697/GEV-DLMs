data {
  int<lower=0> n;
  int<lower=0> k;
  vector [n] y;
  vector[k] FF;
  matrix[k,k] G;
  vector[k] m0;
  matrix[k,k] C0;
}
parameters {
  real xi;
  real<lower=0> tau;
  real<lower=0> scale;
  real<lower=0> sigma;
  vector[n] eta;
}
transformed parameters{
  vector[n] g_t;
  
  matrix[k,k] W; 
 // matrix[d,d] V; 
  
  // DLM Transcitions
  vector[k] a;
  real error;
    
  //array of vectors for the mean
  array[n+1] vector [k] mt;
  vector [n] ft;
  
  vector [k] A;
  matrix[k,k] R;
    
  array[n+1] matrix[k,k] Ct;
  vector<lower=0>[n] Qt;
  
  W =  pow(tau,2)*diag_matrix(rep_vector(1,k));
  
  mt[1] = m0;
  Ct[1] = C0;
  
  for(i in 1:n){
    g_t[i] = scale*pow(xi,-1)*(exp(xi*eta[i]) -1);
    
    a = G * mt[i];
    R = (G *Ct[i] *G') + W;
    
    // Marginal likelihood: 
    ft[i] = sum(FF .* a);
    Qt[i] = (FF' * R * FF) + pow(sigma,2);
    
    // additional computations at t
    error = y[i] - ft[i] - g_t[i];
    A = R * FF / Qt[i];
    
    //  Posterior at t
    mt[i+1] = a + (error*A);
    Ct[i+1] = R - (A  * Qt[i] * A');
  }
}
model {
  // priors
  target += normal_lpdf(xi | 0, 1);
  target += normal_lpdf(scale | 0, 1);
  target += normal_lpdf(sigma | 0, 0.1);
  
  // likelihood
  target += gumbel_lpdf(eta | 0, 1);
  for (i in 1:n) {
     target += normal_lpdf(y[i] - g_t[i] | ft[i], Qt[i]);
   }
}
generated quantities{
  array[n+1] vector[k] theta;
  vector[n] mu;
  vector[n] y_pred;
  matrix[k,k] At;
  matrix[k,k] Q;
  matrix[k,k] H;
  vector[k] h;

  theta[n+1] = multi_normal_rng(mt[n+1],Ct[n+1]);
 
  for(i in 1:n){
    At = Ct[n-i+1] * G';
    Q = inverse(G*At + W);
    H = Ct[n-i+1] - At*Q*At';
    
    h = mt[n-i+1] + At*Q*(theta[n-i+2] - G*mt[n-i+1]);
    theta[n-i+1] = multi_normal_rng(h,H);
    mu[n-i+1] = sum(FF.*theta[n-i+1]);
    y_pred[n-i+1] = mu[n-i+1] + g_t[n-i+1] + normal_rng(0,sigma);
  }
}
