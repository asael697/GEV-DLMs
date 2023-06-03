#' Forward Filtering Backward Sampling algorithm
#' simulates the states parameters in Dynamic linear model
#' 
#'   yt = FF mu_t + et, et ~ N(0,W)
#'   mu_t = G mu_{t-1} + vt, vt ~ N(0,V)
#'   y0 ~ N(m0,C0)
#' The function computes the distribution of the states parameters
#' p(mu_t | y1,y2,...,yt)
#' 
#' @param m0 is the mean of the initial value y0
#' @param C0 is the covariance matrix of the initial value y0
#' @param FF is the link matrix between states (mu_t) and observations (yt)
#' @param G the transition matrix in the states equation.
#' @param n is the length of the data set n = length(y_t)
#' @param m is the dimension of the data yt in R^m
#' @param p is the dimension of the state parameters mu_t in R^p.
#' @param data a matrix of dimensions m rows and n cols containing the data set.
#' @param samples an integer with the amount of samples for approximating the
#' latent distribution of mu_t | y_t.
#' 
#' @author Elvis Arrazola y Cristian Cruz
#' 
#' @return an array mu of dimension n X samples x p
#' 
FFBS <- function(m0 = 0, C0 = 0.6, FF = 1, G = 1, V = 1, 
                 W = 1, n, m = 1, p = 1, data, samples){
  
  at <- mt <- mts <- Mt <-  theta_t <- matrix(ncol = n,nrow = p)
  Rt <- Ct <- Cts <- Bts <- array(dim = c(n,p,p))
  ft <- et <-h_t <- H_t <- matrix(ncol = n,nrow = m)
  Qt <- matrix(ncol = n,nrow = m)
  At <- array(dim = c(n,p,n))
  Theta_t_mayus <- theta_t <- matrix(0,ncol = n,nrow = p)
  
  ##	ECUACIONES ACTUALIZACION PARA t = 1, ... ,length(Go).
  at[,1] = G%*%m0
  Rt[1,,] = (G%*%C0%*%t(G))+W
  ft[,1] = FF%*%at[,1]
  Qt[,1] = FF%*%Rt[1,,]%*%t(FF) + V
  At[1,,] = Rt[1,,]%*%t(FF)%*%solve(Qt[,1])
  et[,1] = data[,1] - ft[,1]  
  mt[,1] = at[,1] + At[1,,]%*%t(et[,1])
  Ct[1,,] = Rt[1,,] - At[1,,]%*%t(Qt[1,1])%*%t(At[1,,])
  
  
  #forward Filter
  for(t in 2:n){
    at[,t] = G%*%mt[,t-1]    
    Rt[t,,] = (G%*%Ct[t-1,,]%*%t(G)) + W
    ft[,t] = FF%*%at[,t]
    Qt[,t] = FF%*%Rt[t,,]%*%t(FF) + V
    At[t,,] = Rt[t,,]%*%t(FF)%*%solve(Qt[,t])
    et[,t] = data[,t] - ft[,t]  
    mt[,t] = at[,t] + At[t,,]%*%t(et[,t])
    Ct[t,,] = Rt[t,,] - At[t,,]%*%t(Qt[,t])%*%t(At[t,,])
    
    
    theta_t[,t] = rnorm(1 , mean =at[,t], sd = sqrt(Rt[t,,]))
    M_t[,t] = rnorm(1, mean = ft[,t] , sd = sqrt(Qt[,t]))
    Theta_t_mayus[,t] = rnorm(1 , mean = M_t[,t] , sd = sqrt(Ct[t,,]))
    
    H_t[,t] = solve(solve(Ct[t,,])+t(G)%*%solve(W)%*%G)
    h_t[,t] = H_t[,t]%*%(solve(Ct[t,,])%*%mt[,t] + t(G)%*%solve(W)%*%M_t[,t])
  }
  # BACKWARD SAMPLING
  if(samples == 1){
    mu = rep(0,n)
    mu[n] = rnorm(1,at[,n],sqrt(Ct[n,,]))
    for (t in (n-1):1) {
      mu[t] = rnorm(1,h_t[,(t+1)],sqrt(H_t[,(t+1)]))
    }
  }else{
    #mus = matrix(0,samples,T)
    mu <- array(dim = c(n,samples,p))
    mu[n,,] = rnorm(samples,at[,n],sqrt(Ct[n,,]))
    
    for (t in (n-1):1) {
      mu[t,,] = rnorm(samples,h_t[,(t+1)],sqrt(H_t[,(t+1)]))
    }
  }

  return(mu)
} 

#' Metropolis Sampler algorithm
#' 
#' Generates different MARKOV Chains with the simulated marginal posteriors
#' 
#' @param y a matrx with the used data for estimating the posterior
#' @param chains integer with the total amount of chains 
#' @param iter an integer with the total iterations per chain after warm up
#' @param scale scale matrix for the proposed Jump distribution
#' @param burnin an integer for the first iterations to be burned in warm-up
#' iterations.
#' @param lag the amount of lags to be burnned for avoid stucked jumps
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return an data frame with the simulated posteriors
#' 
metropolis_sampler <- function(y,chains = 4, iter = 5000, scale = 0.5,
                               burnin = 5000, lag = 5) {
  
  post = NULL
  for(k in 1:chains){
    results = NULL
    current_state = inits()
    # burn-in
    for(i in 1:burnin) {
      out = metropolis_step(y,current_state,scale)
      current_state = out$value
    }
    # sample
    for(i in 1:iter) {
      cont = 0
      for(j in 1:lag) {
        out = metropolis_step(y,current_state,scale)
        current_state = out$value
        cont = cont + out$reject
      }
      results = rbind(results,c(out$value, out$reject))
    }
   post = rbind(post,results) 
  }
  row.names(post) = NULL
  post = as.data.frame(post)
  post$.chain = sort(rep(1:chains,iter))
  
  return(post)
}

#' Metropolis Step
#' 
#' Generates a Jump for the Metropolis Algorithm 
#' 
#' @param y a matrx with the used data for estimating the posterior.
#' @param prop a vector with the previous step
#' @param scale scale matrix for the proposed Jump distribution
#' 
#' This function use a symmetric multivariate normal distribution
#' to generates the Jumps, in case of an asymmetric jump, use the
#' the Metropolis-Hastings jump.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector with the proposed jump
#' 
metropolis_step <- function(y,prop, scale) {
  reject = 1
  proposed = rjump(prop = prop,scale = scale)

  u  =  runif(1)
  t1 = target(y,proposed) + djump(proposed,prop,scale) 
  t2 = target(y,prop) + djump(prop,proposed, scale) 
  et = exp(t1 - t2)
  et = ifelse(is.na(et),0,et)
  accept_prob = min(1,et)
  
  if(u <= accept_prob){
    value = proposed
    reject = 0
  }else{
    value = prop
  }
  
  out = list(value = value,reject = reject)
  return(out)
}

#' Evaluates the target function in a metropolis step
#' 
#'  target(y,theta) = log_lik(y,theta) + log_prior(theta)
#'  
#' 
#' @param y a matrx with the used data for estimating the posterior.
#' @param theta a vector of the unknown parameters to be evaluated in 
#' the prior density.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a real value with the image of the target distribution
#' 
target <- function(y,theta){
  
  t = loglik(y,theta) + log_prior(theta)
  return(t)
}

#' Computes the Log_likelihood matrix 
#' 
#' @param y a matrx with the used data for estimating the posterior.
#' @param post a matrix with the posterior parameters, with ncols = d,
#' and nrows = iter.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a matrix of dimension ncols n = length(y) and nrows = iter,
#' with the evaluated log_likelihood.
#' 
log_lik <- function(y,post){
  
  LL = apply(post, 1, function(z){
    z = as.numeric(z)
    gev_lpdf(y,mu = z[1],sigma = z[2],k = z[3])
  })
  t(LL)
}
#' Defaults log_prior density
#' 
#' @param theta a matrix with the posterior parameters, with ncols = d,
#' and nrows = iter.
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a real value with the image of the log_prior density
#' 
log_prior <- function(theta){
  return(0)
}

#' Compute the density for a GEV distribution
#' 
#' @param y a vector containing the data.
#' @param mu a real value containing the location parameter
#' @param sigma a positive value containing the scale parameter
#' @param k a real value different from zero containing the shape
#' parameter
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a real value with the image of the logarithm of the likelihood
#' 
gev_pdf <-function(y,mu,sigma,k){
  n = length(y)
  z = (y - mu)/sigma
  d = rep(0,n)
  
  if(k == 0){
    d = exp(-z)*exp(-exp(-z))
  }else{
    inv_k = 1/k
    z1 = (1 + k*z)
    cond = k*z > -1
    
    d[cond] = (z1[cond])^(-1-inv_k)*exp(-z1[cond]^(-inv_k))
    d[!cond] = 0
  }
  return(d)
}

#' Simulate values from a GEV distribution
#' 
#' @param y a vector containing the data.
#' @param mu a real value containing the location parameter
#' @param sigma a positive value containing the scale parameter
#' @param k a real value different from zero containing the shape
#' parameter
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector containing a sample that follows a GEV distribution
#' 
gev_rng <- function(n,mu = 0,sigma = 1,k = 1){

  u = runif(n)
  if(k == 0){
    x = mu - sigma*log(-log(u))
  }else{
    x = mu + (sigma/k)*( (-log(u))^k -1)
  }
  return(x)
}

#' Compute the log density for a GEV distribution
#' 
#' @param y a vector containing the data.
#' @param mu a real value containing the location parameter
#' @param sigma a positive value containing the scale parameter
#' @param k a real value different from zero containing the shape
#' parameter
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a real value with the image of the logarithm of the likelihood
#'
gev_lpdf <- function(y,mu,sigma,k){
  n = length(y)
  z = (y - mu)/sigma
  d = rep(0,n)
  
  if(k == 0){
    d = -z - exp(-z)
  }else{
    inv_k = 1/k
    z1 = k*z
    cond = z1 > -1
    
    d[cond] = (-1-inv_k)*log1p(z1[cond]) - (1+z1[cond])^(-inv_k)
    d[!cond] = -1e-64
  }
  return(d)
}

#' Generate a proposed value from a Jump distribution
#' 
#' @param prop a vector with the proposed parameters from previous
#' iteration.
#' @param scale a covariance matrix used in the Gaussian Jump distribution
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector following a Gaussian distribution of dimension d,
#' with mean = prop, and covariance matrix C = scale.
#' 
rjump <- function(prop,scale){
  
  d = length(prop)
  cov = chol(scale)

  x = as.numeric(cov %*%rnorm(n = d))
  
  return(x + prop)
}

#' Evaluate the density of a Gaussian distribution used as Jump density
#' 
#' @param proposed a vector to evaluate the multivariate density.
#' @param prop a vector with the proposed parameters from previous
#' iteration.
#' @param scale a covariance matrix used in the Gaussian Jump distribution
#' 
#' @author Asael Alonzo Matamoros
#' 
#' @return a vector following a Gaussian distribution of dimension d,
#' with mean = prop, and covariance matrix C = scale.
#' 
djump <- function(proposed,prop,scale){
  d2 = mvtnorm::dmvnorm(proposed,mean = prop,
                        sigma = scale,log = TRUE)
  
  return(sum(d2))
}
